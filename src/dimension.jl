
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

# dimension

"""
    Dimension{D} <: AbstractModule

Physical `Dimension` represented by `Group` element `D`.
```Julia
F, M, L, T, Q, Θ, N, J, A, Λ, C
```
Derived `Dimension` can be obtained from multiplicative base.
"""
struct Dimension{D} <: AbstractModule
    @pure Dimension{D}() where D = new{D}()
end

@pure Dimension(D::Group) = Dimension{D}()
@pure isunknown(x::Dimension{D}) where D = isunknown(D)
@pure isdimension(x::Dimension{D}) where D = typeof(D)<:Dimension
@pure dimension(::Dimension{D}) where D = D

value(::Dimension{D}) where D = value(D)

showgroup(io::IO,x::Dimension{D},u) where D = showgroup(io,D,u,'𝟙')
showgroup(io::IO,x,u,c) = showgroup(io,x,u)

Base.show(io::IO,x::Dimension{D}) where D  = show(io,D)

Base.:(==)(::Dimension{A},::Dimension{B}) where {A,B} = A == B

preferred(a,b) = a
function Base.:+(a::Dimension{A},b::Dimension{B}) where {A,B}
    if A.v==B.v
        preferred(a,b)
    elseif islog(A) && islog(B)
        Dimension{A+B}()
    else
        error("addition of Group $a + $b is not valid")
    end
end
function Base.:-(a::Dimension{A},b::Dimension{B}) where {A,B}
    islog(A) && islog(B) ? Dimension{A-B}() : a+b
end

Base.sqrt(::Dimension{D}) where D = Dimension{sqrt(D)}()
Base.inv(::Dimension{D}) where D = Dimension{inv(D)}()
Base.log(x::Dimension{D}) where D = Dimension{log(D)}()
Base.log2(x::Dimension{D}) where D = Dimension{log2(D)}()
Base.log10(x::Dimension{D}) where D = Dimension{log10(D)}()
Base.log(b::Number,x::Dimension{D}) where D = Dimension{log(b,D)}()
Base.exp(x::Dimension{D}) where D = Dimension{exp(D)}()
Base.exp2(x::Dimension{D}) where D = Dimension{exp2(D)}()
Base.exp10(x::Dimension{D}) where D = Dimension{exp10(D)}()
Base.:^(b::T,x::Dimension{D}) where {T<:Number,D} = Dimension{b^D}()
Base.:^(::Dimension{D},b::Number) where D = Dimension{D^b}()
Base.:^(::Dimension{D},b::Integer) where D = Dimension{D^b}()
Base.:^(::Dimension{D},b::Rational) where D = Dimension{D^b}()
Base.:*(::Dimension{A},::Dimension{B}) where {A,B} = Dimension{A*B}()
Base.:/(::Dimension{A},::Dimension{B}) where {A,B} = Dimension{A/B}()
Base.:*(a::Number,b::Dimension{D}) where D = Dimension{a*D}()
Base.:*(a::Dimension{D},b::Number) where D = Dimension{D*b}()
Base.:/(a::Number,b::Dimension) = a*inv(b)
Base.:/(a::Dimension,b::Number) = a*inv(b)

const isq = Values('F','M','L','T','Q','Θ','N','J','A','Λ','C')
const dims = length(isq)
factorfind(x,k,i=0) = iszero(x) ? (x,0) : (r = x%k; iszero(r) ? factorfind(x÷k,k,i+1) : (x,i))

const 𝟙 = Dimension(valueat(0,dims))
for i ∈ 1:dims
    @eval const $(Symbol(isq[i])) = Dimension(valueat($i,dims))
end
const usq = Values(F,M,L,T,Q,Θ,N,J,A,Λ,C)

@pure function factorize(x::Float64,A::Group{Int,N},B::Group{Int,N},C::Group{Int,N},D::Group{Int,N},E::Group{Int,N}) where N
    if isinteger(x)
        try
            return factorize(Int(x),A,B,C,D,E)
        catch
        end
    end
    (x,a) = factorfind(x,2π)
    (x,b) = factorfind(x,π/200)
    (x,c) = factorfind(x,π/180)
    Group(zeros(Values{N,Int}),x)*A^a*(A/B^4/D^2)^b*(A/B^3/C^2/D)^c
end

@pure function factorize(x::Int,A::Group{Int,N},B::Group{Int,N},C::Group{Int,N},D::Group{Int,N},E::Group{Int,N}) where N
    (x,b) = factorfind(x,2)
    (x,c) = factorfind(x,3)
    (x,d) = factorfind(x,5)
    (x,e) = factorfind(x,11)
    Group(zeros(Values{N,Int}),x)*B^b*C^c*D^d*E^e
end

const basis = Values("kB", "NA", "𝘩", "𝘤", "𝘦", "Kcd", "ΔνCs", "R∞", "α", "μₑᵤ", "μₚᵤ", "ΩΛ", "H0", "g₀", "aⱼ", "au", "ft", "ftUS", "lb", "T₀", "atm", "inHg", "RK90", "KJ90", "RK", "KJ", "Rᵤ2014", "Ωᵢₜ", "Vᵢₜ", "kG", "mP", "GME", "GMJ", "φ", "ℯ", "γ", "τ", "2", "3", "5", "7", "11", "19","43")
const vals = length(basis)
@pure function constant(d::Group,C::Coupling=UnitSystems.Universe); cs = 
    UnitSystems.kB^makeint(d.v[1])*
    UnitSystems.NA^makeint(d.v[2])*
    UnitSystems.𝘩^makeint(d.v[3])*
    UnitSystems.𝘤^makeint(d.v[4])*
    UnitSystems.𝘦^makeint(d.v[5])*
    UnitSystems.Kcd^makeint(d.v[6])*
    UnitSystems.ΔνCs^makeint(d.v[7])*
    UnitSystems.R∞^makeint(d.v[8])*
    inv(UnitSystems.αinv)^makeint(d.v[9])*
    UnitSystems.μₑᵤ^makeint(d.v[10])*
    UnitSystems.μₚᵤ^makeint(d.v[11])*
    UnitSystems.ΩΛ^makeint(d.v[12])*
    UnitSystems.H0^makeint(d.v[13])*
    UnitSystems.g₀^makeint(d.v[14])*
    UnitSystems.aⱼ^makeint(d.v[15])*
    UnitSystems.au^makeint(d.v[16])*
    UnitSystems.ft^makeint(d.v[17])*
    UnitSystems.ftUS^makeint(d.v[18])*
    UnitSystems.lb^makeint(d.v[19])*
    UnitSystems.T₀^makeint(d.v[20])*
    UnitSystems.atm^makeint(d.v[21])*
    UnitSystems.inHg^makeint(d.v[22])*
    UnitSystems.RK1990^makeint(d.v[23])*
    UnitSystems.KJ1990^makeint(d.v[24])*
    UnitSystems.RK2014^makeint(d.v[25])*
    UnitSystems.KJ2014^makeint(d.v[26])*
    UnitSystems.Rᵤ2014^makeint(d.v[27])*
    UnitSystems.Ωᵢₜ^makeint(d.v[28])*
    UnitSystems.Vᵢₜ^makeint(d.v[29])*
    UnitSystems.kG^makeint(d.v[30])*
    UnitSystems.mP^makeint(d.v[31])*
    UnitSystems.GME^makeint(d.v[32])*
    UnitSystems.GMJ^makeint(d.v[33])*
    Base.MathConstants.φ^makeint(d.v[34])*
    Base.MathConstants.γ^makeint(d.v[35])*
    Base.MathConstants.ℯ^makeint(d.v[36])*
    (2π)^makeint(d.v[37]); is =
    2.0^makeint(d.v[38])*
    3.0^makeint(d.v[39])*
    5.0^makeint(d.v[40])*
    7.0^makeint(d.v[41])*
    11.0^makeint(d.v[42])*
    19.0^makeint(d.v[43])*
    43.0^makeint(d.v[44])
    return cs*(is*d.c)
end

# convertunit

struct ConvertUnit{D,U,S} <: AbstractModule
    @pure ConvertUnit{D,U,S}() where {D,U,S} = new{D,normal(U),normal(S)}()
end

convertdim(::ConvertUnit{D,U,S}) where {D,U,S} = convertdim(D,U,S)
convertdim(::Dimension{D},U,S) where D = convertdim(D,U,S)
dimconvert(x::T,d,u,s) where T = (isone(ratio(d,u,s)) ? zero(x) : x)::T
convertdim(d::Group{T},U,S) where T = Dimension{Group{T,11}(dimconvert.(d.v,usq,Ref(U),Ref(S)))}()

function Base.show(io::IO,::ConvertUnit{D,U,S}) where {D,U,S}
    d = convertdim(D,U,S)
    print(io, ratio(D,U,S), " [")
    showgroup(io,S(d),S)
    print(io, "]/[")
    showgroup(io,U(d),U)
    print(io, "] ", unitname(U), " -> ", unitname(S))
end

@pure Base.inv(::ConvertUnit{D,U,S}) where {D,U,S} = ConvertUnit{inv(D),U,S}()

(D::Dimension)(U::UnitSystem,S::UnitSystem) = ConvertUnit{D,U,S}()
(D::Dimension)(U::UnitSystem) = U(ratio(D,Natural,U),D)
(D::Dimension)(v::Real,U::UnitSystem,S::UnitSystem=Metric) = v/ratio(D,U,S) #U(_,D)
#(D::Dimension)(v::Number,U::UnitSystem,S::UnitSystem=Natural) = U(v/ratio(D,U,S),D)

Base.log(x::ConvertUnit{D,U,S}) where {D,U,S} = ConvertUnit{log(D),U,S}()
Base.log2(x::ConvertUnit{D,U,S}) where {D,U,S} = ConvertUnit{log2(D),U,S}()
Base.log10(x::ConvertUnit{D,U,S}) where {D,U,S} = ConvertUnit{log10(D),U,S}()
Base.log(b::Number,x::ConvertUnit{D,U,S}) where {D,U,S} = ConvertUnit{log(b,D),U,S}()
Base.exp(x::ConvertUnit{D,U,S}) where {D,U,S} = ConvertUnit{exp(D),U,S}()
Base.exp2(x::ConvertUnit{D,U,S}) where {D,U,S} = ConvertUnit{exp2(D),U,S}()
Base.exp10(x::ConvertUnit{D,U,S}) where {D,U,S} = ConvertUnit{exp10(D),U,S}()
Base.:^(a::Number,b::ConvertUnit{D,U,S}) where {D,U,S} = ConvertUnit{a^D,U,S}()
Base.:^(a::ConvertUnit{D,U,S},b::Integer) where {D,U,S} = Quantity{D^b,U,S}()
Base.:^(a::ConvertUnit{D,U,S},b::Rational{Int}) where {D,U,S} = Quantity{D^b,U,S}()
Base.:*(a::ConvertUnit{A,U,S},b::ConvertUnit{B,U,S}) where {A,B,U,S} = ConvertUnit{A*B,U,S}()
Base.:/(a::ConvertUnit{A,U,S},b::ConvertUnit{B,U,S}) where {A,B,U,S} = ConvertUnit{A/B,U,S}()

# quantity

struct Quantity{D,U,T} <: AbstractModule
    v::T
    Quantity{D,U,T}(v) where {T,D,U} = new{D,normal(U),T}(v)
    @pure Quantity{D,U,Int}(v) where {D,U} = new{D,normal(U),Int}(v)
    @pure Quantity{D,U,Float64}(v) where {D,U} = new{D,normal(U),Float64}(v)
end

Quantity{D,U}(v::Quantity) where {D,U} = Quantity{D,U}(v.v)
Quantity{D,U}(v::T) where {D,U,T} = Quantity{D,U,T}(v)
@pure Quantity{D,U}(v::Int) where {D,U} = Quantity{D,U,Int}(v)::Quantity{D,normal(U),Int}
@pure Quantity{D,U}(v::Float64) where {D,U} = Quantity{D,U,Float64}(v)::Quantity{D,normal(U),Float64}
@pure Quantity(D,U,v::Int) = Quantity{D,U}(v)
@pure Quantity(D,U,v::Float64) = Quantity{D,U}(v)
Quantity(D,U,v) = Quantity{D,U}(v)

Base.convert(::Type{Float64},x::Quantity) = convert(Float64,x.v)
Base.convert(::Type{Float64},x::Quantity{D,U,Float64}) where {D,U} = x.v
Base.convert(::Type{T},x::Quantity{D,U,T}) where {D,U,T<:Number} = x.v

import UnitSystems: evaldim
UnitSystems.isquantity(U::Quantity) = true

quantity(q) = q
quantity(q::Quantity) = q.v
@pure quantity(q::Quantity{D,U,Int} where {D,U}) = q.v
@pure quantity(q::Quantity{D,U,Float64} where {D,U}) = q.v
@pure dimensions(::Quantity{D}) where D = D
@pure Dimension(::Quantity{D}) where D = D
@pure unitsystem(::Quantity{D,U}) where {D,U} = dimension(U)
@pure unitsystem2(::Quantity{D,U}) where {D,U} = U

function Base.show(io::IO,x::Quantity{D,U}) where {D,U}
    print(io, x.v, " [")
    showgroup(io,U(D),U)
    print(io, "] ", unitname(U))
end

Base.log(x::Quantity{D,U}) where {D,U} = Quantity{log(D),U}(log(x.v))
Base.log2(x::Quantity{D,U}) where {D,U} = Quantity{log2(D),U}(log2(x.v))
Base.log10(x::Quantity{D,U}) where {D,U} = Quantity{log10(D),U}(log10(x.v))
Base.log(b::Number,x::Quantity{D,U}) where {D,U} = Quantity{log(b,D),U}(log(b,x.v))
Base.exp(x::Quantity{D,U}) where {D,U} = Quantity{exp(D),U}(exp(x.v))
Base.exp2(x::Quantity{D,U}) where {D,U} = Quantity{exp2(D),U}(exp2(x.v))
Base.exp10(x::Quantity{D,U}) where {D,U} = Quantity{exp10(D),U}(exp10(x.v))
Base.:^(a::Number,b::Quantity{D,U}) where {D,U} = Quantity{a^D,U}(a^b.v)
Base.:^(a::Quantity{D,U},b::Number) where {D,U} = Quantity{D^b,U}(a.v^b)
Base.:^(a::Quantity{D,U},b::Integer) where {D,U} = Quantity{D^b,U}(a.v^b)
Base.:^(a::Quantity{D,U},b::Rational{Int}) where {D,U} = Quantity{D^b,U}(a.v^b)
Base.:+(a::Number,b::Quantity{D,U}) where {D,U} = U(D)==𝟙 ? Quantity{D,U}(a+b.v) : throw(error("$(U(D)) ≠ 𝟙 "))
Base.:+(a::Quantity{D,U},b::Number) where {D,U} = U(D)==𝟙 ? Quantity{D,U}(a.v+b) : throw(error("$(U(D)) ≠ 𝟙 "))
Base.:-(a::Number,b::Quantity{D,U}) where {D,U} = U(D)==𝟙 ? Quantity{D,U}(a-b.v) : throw(error("$(U(D)) ≠ 𝟙 "))
Base.:-(a::Quantity{D,U},b::Number) where {D,U} = U(D)==𝟙 ? Quantity{D,U}(a.v-b) : throw(error("$(U(D)) ≠ 𝟙 "))
Base.:-(a::Quantity{D,U},b::Quantity{D,U}) where {D,U} = Quantity{D,U}(a.v-b.v)
Base.:+(a::Quantity{D,U},b::Quantity{D,U}) where {D,U} = Quantity{D,U}(a.v+b.v)
Base.:*(a::Real,b::Quantity{D,U}) where {D,U} = Quantity{D,U}(a*b.v)
Base.:*(a::Quantity{D,U},b::Real) where {D,U} = Quantity{D,U}(a.v*b)
Base.:*(a::Complex,b::Quantity{D,U}) where {D,U} = Quantity{D,U}(a*b.v)
Base.:*(a::Quantity{D,U},b::Complex) where {D,U} = Quantity{D,U}(a.v*b)
Base.:*(a::Quantity{A,U},b::Quantity{B,U}) where {A,B,U} = Quantity{A*B,U}(a.v*b.v)
Base.:/(a::Quantity{A,U},b::Quantity{B,U}) where {A,B,U} = Quantity{A/B,U}(a.v/b.v)
Base.:/(a::Quantity{D,U},b::Quantity{D,U}) where {D,U} = Quantity{D/D,U}(a.v/b.v)
Base.:/(a::Quantity{D,A},b::Quantity{D,B}) where {A,B,D} = ConvertUnit{(a.v/b.v)*D,A,B}()
Base.:/(a::Number,b::Quantity{D}) where D = a*inv(b)
Base.:/(a::Quantity{D},b::Number) where D = a*inv(b)
Base.:/(a::Quantity{D},b::UnitSystems.Constant) where D = a*inv(b)
Base.:/(a::UnitSystems.Constant,b::Quantity{D}) where D = a*inv(b)
Base.:/(a::Quantity{D,U},b::ConvertUnit{D,S,U}) where {D,U,S} = a*ConvertUnit{D,U,S}()
Base.:-(a::Quantity{D,U}) where {D,U} = Quantity{D,U}(-a.v)
Base.inv(a::Quantity{D,U}) where {D,U} = Quantity{inv(D),U}(inv(a.v))
Base.sqrt(a::Quantity{D,U}) where {D,U} = Quantity{sqrt(D),U}(sqrt(a.v))

Base.:*(a::Dimension{A},b::Quantity{D,U}) where {D,U,A} = Quantity{a,U}(10^-(isunknown(A) ? A.v.v[end] : A.v[end]))*b
Base.:*(a::Quantity{D,U},b::Dimension{B}) where {D,U,B} = a*Quantity{b,U}(10^-(isunknown(B) ? B.v.v[end] : B.v[end]))
Base.:*(a::Quantity{D,U},b::ConvertUnit{A,U,U}) where {D,U,A} = D==A || A==𝟙 ? a : error("ConvertUnit incompatible $D ≠ $A")
Base.:*(a::ConvertUnit{A,U,U},b::Quantity{D,U}) where {D,U,A} = D==A || A==𝟙 ? b : error("ConvertUnit incompatible $D ≠ $A")
Base.:*(a::ConvertUnit{A,S,U},b::Quantity{D,U}) where {D,U,A,S} = (A==D && S==U) ? b : inv(A)==D ? ConvertUnit{inv(A),U,S}()*b : error("ConvertUnit incompatible $(inv(A)) ≠ $D")
Base.:*(a::Quantity{D,U},b::ConvertUnit{A,S,U}) where {D,U,A,S} = (D==A && S==U) ? a : inv(A)==D ? a*ConvertUnit{inv(A),U,S}() : error("ConvertUnit incompatible $D ≠ $(inv(A))")
Base.:*(a::ConvertUnit{A,U,S},b::Quantity{D,U}) where {A,D,U,S} = A==D ? Quantity{D,S}(ratio(D,U,S)*b.v) : error("ConvertUnit incompatible $A ≠ $B")
Base.:*(a::Quantity{D,U},b::ConvertUnit{A,U,S}) where {A,D,U,S} = A==D ? Quantity{D,S}(a.v*ratio(D,U,S)) : error("ConvertUnit incompatible $A ≠ $B")

Base.:(==)(a::Quantity{A,U},b::Quantity{B,U}) where {A,B,U} = U(A) == U(B) && a.v == b.v
Base.:(==)(a::Number,b::Quantity{D}) where D = iszero(norm(first(value(D),dims-1))) && a == b.v*10^last(value(D))
Base.:(==)(a::Quantity{D},b::Number) where D = iszero(norm(first(value(D),dims-1))) && b == a.v*10^last(value(D))

@pure isunknown(x) = false

@pure sumabs(::Dimension{D}) where D = sum(abs.(D.v))
@pure add(d::Dimension{D},::Dimension{D},U) where D = islog(D) ? d+d : d
@pure function add(a::Dimension{A},b::Dimension{B},U) where {A,B}
    Ua,Ub = U(a),U(b)
    if Ua == Ub
        sa,sb = sumabs(a),sumabs(b)
        sb < sa || (sb == sa && Ub == b) ? b : a
    elseif islog(A) && islog(B)
        a + b
    else
        throw(error("$(Ua) ≠ $(Ub)"))
    end
end
@pure sub(d::Dimension{D},::Dimension{D},U) where D = islog(D) ? 𝟙 : d
@pure function sub(a::Dimension{A},b::Dimension{B},U) where {A,B}
    islog(A) && islog(B) ? a - b : add(a,b)
end

Base.:+(a::Quantity{A,U},b::Quantity{B,U}) where {A,B,U} = Quantity{add(A,B,U),U}(a.v+b.v)
Base.:-(a::Quantity{A,U},b::Quantity{B,U}) where {A,B,U} = Quantity{sub(A,B,U),U}(a.v-b.v)

# other

Base.isone(x::Dimension{D}) where D = isone(D)
Base.isone(x::Quantity{D}) where D = false
Base.iszero(x::Dimension{D}) where D = iszero(D)

# Quantities

struct Quantities{D,U,N,T} <: TupleVector{N,T}
    v::Values{N,T}
    Quantities{D,U,N,T}(v::Values{N,T}) where {D,U,N,T} = new{D,normal(U),N,T}(v)
end

Quantities{D,U,N,T}(v::Tuple) where {D,U,N,T} = Quantities{D,U,N,T}(Values{N,T}(v))
Quantities{D,U}(v::Values{N,T}) where {D,U,N,T} = Quantities{D,U,N,T}(v)
Quantities(D,U,v) = Quantities{D,U}(v)

(S::UnitSystem)(x::Quantities{D,U}) where {D,U} = Quantities{D,S}(D.(x.v,Ref(U),Ref(S)))

export Quantity, Quantities

Base.getindex(x::Quantities{D,U},i::Integer) where {D,U} = Quantity{D,U}(x.v[i])
Base.getindex(x::Quantities{D,U},i::Int) where {D,U} = Quantity{D,U}(x.v[i])

# functors

import UnitSystems: UnitSystem
UnitSystem(::Dimension{D}) where D = Dimension{UnitSystem(D)}()
UnitSystem(g::LogGroup{B}) where B = LogGroup{B}(UnitSystem(value(g)))
UnitSystem(g::ExpGroup{B}) where B = ExpGroup{B}(UnitSystem(value(g)))
function UnitSystem(d::Group)
    Group(Values(
        -d.v[6],
        d.v[3]+d.v[4]+d.v[5]/2-d.v[1]-d.v[8],
        3(d.v[1])+2d.v[6]+4d.v[8]-d.v[3]-2d.v[4]-d.v[5]/2,
        d.v[5]/-2,
        d.v[2]+d.v[6]+d.v[7]+2(d.v[1]+d.v[8])-d.v[3]-d.v[4],
        -d.v[7],
        d.v[8],
        d.v[9]-2d.v[10],
        d.v[10]-d.v[5]/2,
        -(d.v[5]+d.v[11]),
        d.v[3]+d.v[4]-d.v[6]-2(d.v[1]+d.v[8])))
end
function UnitSystem(d::Group{<:Integer})
    Group(Values(
        -d.v[6],
        d.v[3]+d.v[4]+d.v[5]//2-d.v[1]-d.v[8],
        3(d.v[1])+2d.v[6]+4d.v[8]-d.v[3]-2d.v[4]-d.v[5]//2,
        d.v[5]//-2,
        d.v[2]+d.v[6]+d.v[7]+2(d.v[1]+d.v[8])-d.v[3]-d.v[4],
        -d.v[7],
        d.v[8],
        d.v[9]-2d.v[10],
        d.v[10]-d.v[5]//2,
        -(d.v[5]+d.v[11]),
        d.v[3]+d.v[4]-d.v[6]-2(d.v[1]+d.v[8])))
end

(u::UnitSystem)(a::Number, d::Dimension) = Quantity{d,u}(a)
(u::UnitSystem)(a::Number, d::AbelianGroup) = Quantity{Dimension(d),u}(a)
(u::UnitSystem)(::Dimension{D}) where D = Dimension{normal(u)(D)}()
(u::UnitSystem)(d::LogGroup{B}) where B = LogGroup{B}(u(d.v))
(u::UnitSystem)(d::ExpGroup{B}) where B = ExpGroup{B}(u(d.v))
