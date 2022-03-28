
#   This file is part of Similitude.jl
#   It is licensed under the MIT license
#   Similitude Copyright (C) 2022 Michael Reed
#       _           _                         _
#      | |         | |                       | |
#   ___| |__   __ _| | ___ __ __ ___   ____ _| | __ _
#  / __| '_ \ / _` | |/ / '__/ _` \ \ / / _` | |/ _` |
# | (__| | | | (_| |   <| | | (_| |\ V / (_| | | (_| |
#  \___|_| |_|\__,_|_|\_\_|  \__,_| \_/ \__,_|_|\__,_|
#
#   https://github.com/chakravala
#   https://crucialflow.com

abstract type AbstractModule <: Number end

abstract type AbelianGroup <: AbstractModule end

struct Group{T,N,S} <: AbelianGroup
    v::Values{N,T}
    c::S
    @pure Group{T,N,S}(v,c) where {N,T,S} = new{T,N,S}(v,UnitSystems.cache(c))
end

@pure Group{T,N}(v,c::S=1) where {T,N,S} = Group{T,N,S}(v,c)
@pure Group{T}(v::Values{N,T},c::S=1) where {N,T,S} = Group{T,N,S}(v,c)
@pure Group(v::Values{N,T},c::S=1) where {N,T,S} = Group{T,N,S}(v,c)

@pure function Group(v::Values{N,<:Rational},c::S=1) where {N,S}
    r = promoteint(v)
    Group{eltype(r),N,S}(r,c)
end

import AbstractTensors: value
value(g::Group) = g.v
coef(g::Group) = UnitSystems.measure(g.c)

Base.:(==)(a::Group,b::Group) = a.v == b.v && a.c == b.c

@pure isonezero(x) = isone(x) || iszero(x)
@pure checkint(v::Values{N} where N) = prod(isonezero.(v.v))
@pure checkint(v::Values{N,<:Rational} where N) = prod(isone.(denominator.(v.v)))
@pure checkint(v::Values{N,<:Integer} where N) = v
@pure promoteint(v::Values{N,<:Integer} where N) = v
@pure promoteint(v) = checkint(v) ? Int.(v) : v

const expos = Values('⁰','¹','²','³','⁴','⁵','⁶','⁷','⁸','⁹')
const chars = Dict([[string(i-1)[1]=>expos[i] for i ∈ 1:length(expos)];['.'=>'⋅','-'=>'⁻','e'=>'ᵉ','v'=>'ᵛ','₀'=>'⁰','₁'=>'¹','₂'=>'²','₃'=>'³','₄'=>'⁴','₅'=>'⁵','₆'=>'⁶','₇'=>'⁷','₈'=>'⁸','₉'=>'⁹','*'=>'*']])

makeint(x) = x
@pure makint(x::Int) = x
@pure function makeint(x::AbstractFloat)
    ax,rem,ne = abs(x),abs(x%1),sqrt(eps()*norm(x))
    if ne < 1
        if log10(ax)-log(1.7,rem) > 20
            return Int(x÷1)
        elseif log10(ax)-log10(1-rem) > 17
            return Int(x÷1)+1
        else
            return x
        end
    else
        return x
    end
end

findpower(x) = 0
findpower(x::Int,i=Int(round(log10(x),RoundDown))) = i<0 ? 0 : (d = (x÷(10^i))%10; !iszero(d) ? findpower(x,i-1) : i+1)

printexpo(io::IO, d, x) = !iszero(x) && (print(io, d); printexpo(io,x))
function printexpo(io::IO, d, x::AbstractFloat)
    if !iszero(x)
        ix = makeint(x)
        if abs(x) < 1
            mix = makeint(inv(x))
            if typeof(mix) == Int
                printexpo(io, d, 1//mix)
            else
                if (d == 10 || d == "10") && length(string(abs(x)))>5
                    x < 0 && print(io, '/')
                    print(io, makeint(10^abs(x)))
                    !(x<0) && print(io, '*')
                else
                    print(io, d)
                    printexpo(io, x)
                end
            end
        elseif d == 10 || d == "10"
            mx = makeint(abs(x))
            x < 0 && print(io, '/')
            if typeof(mx)== Int
                printexpo(io, d, mx)
            else
                ten = makeint(10^abs(x))
                pow = findpower(ten)
                if !iszero(pow)
                    net = ten÷10^pow
                    if !isone(net)
                        print(io, net)
                        print(io, x < 0 ? '/' : '*')
                    end
                    print(io, d)
                    printexpo(io, pow)
                elseif length(string(abs(x)))>5
                    print(io, ten)
                    !(x<0) && print(io, '*')
                else
                    printexpo(io, d, rationalize(x))
                end
                if !iszero(pow)
                    
                end
            end
        elseif typeof(ix) == Int
            printexpo(io,d,ix)
        else
            print(io, d)
            printexpo(io,x)
        end
    end
end

function printexpo(io::IO, x::Integer)
    if !isone(x)
        print(io, (x<0 ? ('⁻',) : ())..., expos[reverse(digits(abs(x)).+1)]...)
    end
end

function printexpo(io::IO, x::AbstractFloat)
    if !isone(x)
        print(io, (x<0 ? ('⁻',) : ())..., [chars[i] for i ∈ string(abs(x))]...)
    end
end

function printexpo(io::IO, x::Rational)
    if !isone(x)
        print(io, (x.num<0 ? ('⁻',) : ())..., expos[reverse(digits(abs(x.num)).+1)]...)
        !isone(x.den) && print(io, 'ᐟ', expos[reverse(digits(x.den).+1)]...)
    end
end

function printexpo(io::IO, x::Complex)
    if !isone(x)
        if !iszero(x.re)
            isone(x.re) ? print(io, '¹') : printexpo(io, x.re)
            !iszero(x.im) && print(io, x.im<0 ? '⁻' : '⁺')
        else
            x.im<0 && print(io, '⁻')
        end
        if !iszero(x.im)
            printexpo(io, abs(x.im))
            print(io, "ⁱᵐ")
        end
    end
end

function printdims(io::IO,x::Group{T,N},name) where {T,N}
    M = 0
    if haskey(ENV,"GROUPAREN")
        for i ∈ 1:N-M
            printnum(io, name[i], x.v[i])
        end
        n = sum(first(x.v,N-M).<0)
        n>0 && print(io, '/')
        n>1 && print(io, '(')
        for i ∈ 1:N-M
            printnum(io, name[i], -x.v[i])
            typeof(name[i])==String && isone(x.v[i]) && !iszero(norm(last(first(x.v,N-M),N-i-1))) && print(io,'*')
        end
        n>1 && print(io, ')')
    else
        for i ∈ 1:N-M
            printexpo(io, name[i], makeint(x.v[i]))
            typeof(name[i])==String && isone(x.v[i]) && !iszero(norm(last(first(x.v,N-M),N-i-M))) && print(io,'*')
        end
    end
end

function printnum(io, b, e)
    me = makeint(e)
    !(me<0) && printexpo(io, b, me)
end

Base.show(io::IO,x::Group) = showgroup(io,x)
function showgroup(io::IO,x::Group{T,N},c='𝟙',u=Natural) where {T,N}
    #back = T<:AbstractFloat && x.v[N]<0
    #!back && printexpo(io, 10, x.v[N])
    printdims(io,x,u)
    iz = iszero(norm(x.v))
    xc = coef(x)
    iz && (isone(xc)||abs(xc)<1) && print(io, c)
    #back && printexpo(io, 10, last(x.v))
    if !isone(xc)
        if abs(xc)<1
            print(io,'/',makeint(inv(xc)))
        else
            !iz && print(io, '*')
            print(io, makeint(xc))
        end
    end
end

function showgroup(io::IO, x::AbelianGroup, c='𝟙',u=Metric)
    showfun(io,x)
    showgroup(io,dimensions(x),c,u)
    print(io,')')
end

# log

struct LogGroup{B,T<:AbelianGroup} <: AbelianGroup
    v::T
    @pure LogGroup{B,T}(v) where {B,T} = new{B,T}(v)
end

@pure LogGroup{B}(v::T) where {B,T} = LogGroup{B,T}(v)
@pure LogGroup(d::AbelianGroup) = LogGroup{ℯ}(d)

@pure islog(::LogGroup) = true
@pure islog(x) = false

@pure base(x::LogGroup{B}) where B = B
@pure dimensions(x::LogGroup) = x.v

value(g::LogGroup) = g.v

Base.show(io::IO, x::LogGroup) = (showfun(io,x); print(io,dimensions(x),')'))

showfun(io::IO, x::LogGroup{B}) where B = print(io,"log(",B,',')
showfun(io::IO, x::LogGroup{ℯ}) = print(io,"log(")
showfun(io::IO, x::LogGroup{2}) = print(io,"log2(")
showfun(io::IO, x::LogGroup{10}) = print(io,"log10(")
showfun(io::IO, x::LogGroup{exp10(0.1)}) = print(io,"dB(")

Base.log(x::AbelianGroup) = LogGroup(x)
Base.log2(x::AbelianGroup) = LogGroup{2}(x)
Base.log10(x::AbelianGroup) = LogGroup{10}(x)
Base.log(b::Number,x::AbelianGroup) = LogGroup{b}(x)
Base.exp(x::LogGroup{ℯ}) = dimensions(x)
Base.exp2(x::LogGroup{2}) = dimensions(x)
Base.exp10(x::LogGroup{10}) = dimensions(x)
Base.exp(x::LogGroup{B}) where B = dimensions(x)^inv(log(B))
Base.exp2(x::LogGroup) = exp2(x)
Base.exp10(x::LogGroup) = exp10(x)
Base.:^(b::T,x::LogGroup) where T<:Number = exp(x*log(b))

Base.:+(x::LogGroup{B},y::LogGroup{B}) where B = LogGroup{B}(x.v*y.v)
Base.:-(x::LogGroup{B},y::LogGroup{B}) where B = LogGroup{B}(x.v/y.v)
Base.:/(x::LogGroup{B},y::T) where {B,T<:Number} = LogGroup{B^y}(x.v)
Base.:*(x::LogGroup,y::T) where T<:Number = x/inv(y)
Base.:*(x::T,y::LogGroup) where T<:Number = y*x

# exp

struct ExpGroup{B,T<:AbelianGroup} <: AbelianGroup
    v::T
end

@pure ExpGroup{B}(v::T) where {B,T} = ExpGroup{B,T}(v)
@pure ExpGroup(d::AbelianGroup) = ExpGroup{ℯ}(d)

@pure base(x::ExpGroup{B}) where B = B
@pure dimensions(x::ExpGroup) = x.v

value(g::ExpGroup) = g.v

Base.show(io::IO, x::ExpGroup) = (showfun(io,x); print(io,dimensions(x),')'))

showfun(io::IO, x::ExpGroup{B}) where B = print(io,B,"^(")
showfun(io::IO, x::ExpGroup{ℯ}) = print(io,"exp(")
showfun(io::IO, x::ExpGroup{2}) = print(io,"exp2(")
showfun(io::IO, x::ExpGroup{10}) = print(io,"exp10(")

Base.exp(x::AbelianGroup) = ExpGroup(x)
Base.exp2(x::AbelianGroup) = ExpGroup{2}(x)
Base.exp10(x::AbelianGroup) = ExpGroup{10}(x)
Base.:^(b::T,x::AbelianGroup) where T<:Number = ExpGroup{b}(x)
Base.log(x::ExpGroup{ℯ}) = dimensions(x)
Base.log2(x::ExpGroup{2}) = dimensions(x)
Base.log10(x::ExpGroup{10}) = dimensions(x)
Base.log(b::Number,x::ExpGroup{B}) where B = dimensions(x)/log(B,b)
Base.log(x::ExpGroup) = log(ℯ,x)
Base.log2(x::ExpGroup) = log(2,x)
Base.log10(x::ExpGroup) = log(10,x)

Base.:^(x::ExpGroup{B},y::T) where {B,T<:Number} = iszero(y) ? one(x) : x
Base.:^(x::ExpGroup{B},y::T) where {B,T<:Integer} = iszero(y) ? one(x) : x
Base.:^(x::ExpGroup,y::ExpGroup) = ExpGroup{x}(y)
Base.:*(x::ExpGroup{B},y::ExpGroup{B}) where B = ExpGroup{B}(x.v+y.v)
Base.:/(x::ExpGroup{B},y::ExpGroup{B}) where B = ExpGroup{B}(x.v-y.v)

Base.:*(x::ExpGroup{X},y::ExpGroup{Y}) where {X,Y} = ExpGroup{X*Y}(x.v+y.v)

# other

dB(x::AbelianGroup) = log(exp10(0.1),x)

Base.one(::AbelianGroup) = Group(zeros(Values{dims,Int}))
Base.isone(x::Group) = iszero(norm(x.v)) && isone(x.c)
Base.isone(x::AbelianGroup) = false
Base.zero(x::AbelianGroup) = log(one(x))
Base.iszero(x::AbelianGroup) = false
Base.iszero(x::LogGroup) = isone(dimensions(x))

Base.:*(a::Group,b::Group) = Group(a.v+b.v,coef(a)*coef(b))
Base.:/(a::Group,b::Group) = Group(a.v-b.v,coef(a)/coef(b))
Base.:^(a::Group,b::Number) = Group(b*a.v,coef(a)^b)
Base.:^(a::Group,b::Integer) = Group(b*a.v,coef(a)^b)
Base.:^(a::Group,b::Rational) = Group(b*a.v,coef(a)^b)
Base.sqrt(a::Group{T}) where T = Group{T}(a.v/2,sqrt(coef(a)))
Base.sqrt(a::Group{Int}) = Group{Rational{Int}}(a.v//2,sqrt(coef(a)))
Base.inv(a::Group{T}) where T = Group{T}(-a.v,inv(coef(a)))

@pure valueat(::Val{j},::Val{k}) where {j,k} = valueat(j,k)
valueat(j,k,z::T=1) where T = Group(Values{k,T}([i==j ? z : 0 for i ∈ 1:k]))
Base.:*(a::Number,b::Group{T,N}) where {T,N} = Group{T,N}(b.v,a*coef(b))
Base.:*(a::Group{T,N},b::Number) where {T,N} = Group{T,N}(a.v,coef(a)*b)
Base.:/(a::Number,b::Group) = a*inv(b)#
Base.:/(a::Group,b::Number) = a*inv(b)#
