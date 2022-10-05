
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

import UnitSystems: listext, Kinematic, Mechanical, Electromagnetic, Thermodynamic, Molar, Photometric

# dimension

@pure Constant(D::Group) = Constant{D}()
@pure dimension(::Constant{D}) where D = D

value(::Constant{D}) where D = value(D)

showgroup(io::IO,x::Constant{D},u) where D = showgroup(io,D,dimtext(u),'𝟙')
showgroup(io::IO,x,u,c) = showgroup(io,x,u)

preferred(a,b) = a

function Base.:+(a::Group{:USQ},b::Group{:USQ})
    if a.v==b.v && a.c==b.c
        preferred(a,b)
    #elseif islog(A) && islog(B)
    #    Dimension{A+B}()
    else
        error("addition of Group $a + $b is not valid")
    end
end
function Base.:-(a::Group{:USQ},b::Group{:USQ})
    a+b
    #islog(A) && islog(B) ? Dimension{A-B}() : a+b
end

function Base.:+(a::Group{:Constants},b::Group{:Constants})
    if a.v==b.v
        a.c==b.c ? 𝟐*preferred(a,b) : Group(preferred(a,b).v,a.c+b.c,Val(:Constants))
    #elseif islog(A) && islog(B)
    #    Dimension{A+B}()
    else
        product(a)+product(b)
    end
end
function Base.:-(a::Group{:Constants},b::Group{:Constants})
    if a.v==b.v
        a.c==b.c ? 0 : Group(preferred(a,b).v,a.c-b.c,Val(:Constants))
    #elseif islog(A) && islog(B)
    #    Dimension{A+B}()
    else
        product(a)-product(b)
    end
end
Base.:+(a::Real,b::Group{:Constants}) = a+product(b)
Base.:+(a::Group{:Constants},b::Real) = product(a)+b
Base.:-(a::Real,b::Group{:Constants}) = a-product(b)
Base.:-(a::Group{:Constants},b::Real) = product(a)-b

Base.:*(a::Group{:Constants},b::Group{:USQ}) = Group(b.v,a*b.c,Val(:USQ))
Base.:*(a::Group{:USQ},b::Group{:Constants}) = Group(a.v,a.c*b,Val(:USQ))
Base.:/(a::Group{:Constants},b::Group{:USQ}) = a*inv(b)
Base.:/(a::Group{:USQ},b::Group{:Constants}) = a*inv(b)

@group USQ F M L T Q Θ N J A R C

const isq = Values('F','M','L','T','Q','Θ','N','J','A','R','C')
const dims = length(isq)

const 𝟙 = Constant(valueat(0,dims,:USQ))
const USQ = Values(F,M,L,T,Q,Θ,N,J,A,R,C)
const usq = USQ

@group Constants begin
    kB = UnitSystems.kB
    NA = UnitSystems.NA
    𝘩 = UnitSystems.𝘩
    𝘤 = UnitSystems.𝘤
    𝘦 = UnitSystems.𝘦
    Kcd = UnitSystems.Kcd
    ΔνCs = UnitSystems.ΔνCs
    R∞ = UnitSystems.R∞
    α = UnitSystems.α
    μₑᵤ = UnitSystems.μₑᵤ
    μₚᵤ = UnitSystems.μₚᵤ
    ΩΛ = UnitSystems.ΩΛ
    H0 = UnitSystems.H0
    g₀ = UnitSystems.g₀
    aⱼ = UnitSystems.aⱼ
    au = UnitSystems.au
    ft = UnitSystems.ft
    ftUS = UnitSystems.ftUS
    lb = UnitSystems.lb
    T₀ = UnitSystems.T₀
    atm = UnitSystems.atm
    inHg = UnitSystems.inHg
    RK90 = UnitSystems.RK1990
    KJ90 = UnitSystems.KJ1990
    RK = UnitSystems.RK2014
    KJ = UnitSystems.KJ2014
    Rᵤ2014 = UnitSystems.Rᵤ2014
    Ωᵢₜ = UnitSystems.Ωᵢₜ
    Vᵢₜ = UnitSystems.Vᵢₜ
    kG = UnitSystems.kG
    mP = UnitSystems.mP
    GME = UnitSystems.GME
    GMJ = UnitSystems.GMJ
    φ = Base.MathConstants.φ
    γ = Base.MathConstants.γ
    ℯ = Base.MathConstants.ℯ
    τ ≡ 2π
    2 = 2
    3 = 3
    5 = 5
    7 = 7
    11 = 11
    19 = 19
    43 = 43
end

Base.show(io::IO,x::Group{:Constants}) = showgroup(io,x,basis,'𝟏')
@pure promoteint(v::Constant) = isone(v) ? 1 : v

phys(j,k=vals) = Constant(valueat(j,k,:Constants))

const basis = Values("kB", "NA", "𝘩", "𝘤", "𝘦", "Kcd", "ΔνCs", "R∞", "α", "μₑᵤ", "μₚᵤ", "ΩΛ", "H0", "g₀", "aⱼ", "au", "ft", "ftUS", "lb", "T₀", "atm", "inHg", "RK90", "KJ90", "RK", "KJ", "Rᵤ2014", "Ωᵢₜ", "Vᵢₜ", "kG", "mP", "GME", "GMJ", "φ", "γ", "ℯ", "τ", "2", "3", "5", "7", "11", "19","43")
const vals = length(basis)

# convertunit

struct ConvertUnit{D,U,S} <: AbstractModule
    @pure ConvertUnit{D,U,S}() where {D,U,S} = new{D,normal(U),normal(S)}()
end

convertdim(::ConvertUnit{D,U,S}) where {D,U,S} = convertdim(D,U,S)
convertdim(::Constant{D},U,S) where D = convertdim(D,U,S)
dimconvert(x::T,d,u,s) where T = (isone(ratio(d,u,s)) ? zero(x) : x)::T
convertdim(d::Group{:USQ,T},U,S) where T = Constant{Group{:USQ,T}(dimconvert.(d.v,usq,Ref(U),Ref(S)))}()

function Base.show(io::IO,::ConvertUnit{D,U,S}) where {D,U,S}
    d = convertdim(D,U,S)
    print(io, ratio(D,U,S), " [")
    showgroup(io,S(d),S)
    print(io, "]/[")
    showgroup(io,U(d),U)
    print(io, "] ", unitname(U), " -> ", unitname(S))
end

@pure Base.inv(::ConvertUnit{D,U,S}) where {D,U,S} = ConvertUnit{inv(D),U,S}()

(D::Constant)(U::UnitSystem,S::UnitSystem) = ConvertUnit{D,U,S}()
(D::Constant)(U::UnitSystem) = U(ratio(D,Natural,U),D)
(D::Constant)(v::Real,U::UnitSystem,S::UnitSystem=Metric) = v/ratio(D,U,S) #U(_,D)
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
#Base.:^(a::Constant,b::Quantity{D,U}) where {D,U} = Quantity{a^D,U}(a^b.v)
Base.:^(a::UnitSystems.Constant,b::Quantity{D,U}) where {D,U} = Quantity{a^D,U}(a^b.v)
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
Base.cbrt(a::Quantity{D,U}) where {D,U} = Quantity{cbrt(D),U}(cbrt(a.v))

#Base.:*(a::Dimension{A},b::Quantity{D,U}) where {D,U,A} = Quantity{a,U}(10^-(isunknown(A) ? A.v.v[end] : A.v[end]))*b
#Base.:*(a::Quantity{D,U},b::Dimension{B}) where {D,U,B} = a*Quantity{b,U}(10^-(isunknown(B) ? B.v.v[end] : B.v[end]))
Base.:*(a::Quantity{D,U},b::ConvertUnit{A,U,U}) where {D,U,A} = D==A || A==𝟙 ? a : error("ConvertUnit incompatible $D ≠ $A")
Base.:*(a::ConvertUnit{A,U,U},b::Quantity{D,U}) where {D,U,A} = D==A || A==𝟙 ? b : error("ConvertUnit incompatible $D ≠ $A")
Base.:*(a::ConvertUnit{A,S,U},b::Quantity{D,U}) where {D,U,A,S} = (A==D && S==U) ? b : inv(A)==D ? ConvertUnit{inv(A),U,S}()*b : error("ConvertUnit incompatible $(inv(A)) ≠ $D")
Base.:*(a::Quantity{D,U},b::ConvertUnit{A,S,U}) where {D,U,A,S} = (D==A && S==U) ? a : inv(A)==D ? a*ConvertUnit{inv(A),U,S}() : error("ConvertUnit incompatible $D ≠ $(inv(A))")
Base.:*(a::ConvertUnit{A,U,S},b::Quantity{D,U}) where {A,D,U,S} = A==D ? Quantity{D,S}(ratio(D,U,S)*b.v) : error("ConvertUnit incompatible $A ≠ $B")
Base.:*(a::Quantity{D,U},b::ConvertUnit{A,U,S}) where {A,D,U,S} = A==D ? Quantity{D,S}(a.v*ratio(D,U,S)) : error("ConvertUnit incompatible $A ≠ $B")

Base.:+(a::Constant,b::Quantity{D,U}) where {D,U} = U(D)==𝟙 ? Quantity{D,U}(a+b.v) : throw(error("$(U(D)) ≠ 𝟙 "))
Base.:+(a::Quantity{D,U},b::Constant) where {D,U} = U(D)==𝟙 ? Quantity{D,U}(a.v+b) : throw(error("$(U(D)) ≠ 𝟙 "))
Base.:-(a::Constant,b::Quantity{D,U}) where {D,U} = U(D)==𝟙 ? Quantity{D,U}(a-b.v) : throw(error("$(U(D)) ≠ 𝟙 "))
Base.:-(a::Quantity{D,U},b::Constant) where {D,U} = U(D)==𝟙 ? Quantity{D,U}(a.v-b) : throw(error("$(U(D)) ≠ 𝟙 "))
Base.:*(a::Constant,b::Quantity{D,U}) where {D,U} = Quantity{D,U}(a*b.v)
Base.:*(a::Quantity{D,U},b::Constant) where {D,U} = Quantity{D,U}(a.v*b)
#Base.:/(a::Constant,b::Quantity{D,U}) where {D,U} = Quantity{inv(D),U}(a/b.v)
#Base.:/(a::Quantity{D,U},b::Constant) where {D,U} = Quantity{D,U}(a.v/b)

Base.:(==)(a::Quantity{A,U},b::Quantity{B,U}) where {A,B,U} = U(A) == U(B) && a.v == b.v
Base.:(==)(a::Number,b::Quantity{D}) where D = iszero(norm(first(value(D),dims-1))) && a == b.v*10^last(value(D))
Base.:(==)(a::Quantity{D},b::Number) where D = iszero(norm(first(value(D),dims-1))) && b == a.v*10^last(value(D))

@pure isunknown(x) = false

@pure sumabs(::Constant{D}) where D = sum(abs.(D.v))
@pure add(d::Constant{D},::Constant{D},U) where D = islog(D) ? d+d : d
@pure function add(a::Constant{A},b::Constant{B},U) where {A,B}
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
@pure sub(d::Constant{D},::Constant{D},U) where D = islog(D) ? 𝟙 : d
@pure function sub(a::Constant{A},b::Constant{B},U) where {A,B}
    islog(A) && islog(B) ? a - b : add(a,b)
end

Base.:+(a::Quantity{A,U},b::Quantity{B,U}) where {A,B,U} = Quantity{add(A,B,U),U}(a.v+b.v)
Base.:-(a::Quantity{A,U},b::Quantity{B,U}) where {A,B,U} = Quantity{sub(A,B,U),U}(a.v-b.v)

# other

Base.isone(x::Quantity{D}) where D = false

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
# 1,2,3,4, 5, 6, 7,  8,9,10,11
#kB,ħ,𝘤,μ₀,mₑ,Mᵤ,Kcd,A,λ,αL,g₀
# F,M,L,T, Q, Θ, N,  J,A,R, C

import UnitSystems: UnitSystem
UnitSystem(::Constant{D}) where D = Constant{UnitSystem(D)}()
UnitSystem(g::LogGroup{B}) where B = LogGroup{B}(UnitSystem(value(g)))
UnitSystem(g::ExpGroup{B}) where B = ExpGroup{B}(UnitSystem(value(g)))
function UnitSystem(d::Group{:USQ})
    Group(Values(
        -d.v[6],
        d.v[3]+d.v[4]+d.v[5]/2-d.v[1]-d.v[8],
        3(d.v[1])+2d.v[6]+4d.v[8]-d.v[3]-2d.v[4]-d.v[5]/2,
        d.v[5]/-2,
        d.v[2]+d.v[6]+d.v[7]+2(d.v[1]+d.v[8])-d.v[3]-d.v[4],
        -d.v[7],
        d.v[8],
        d.v[3]+d.v[4]+d.v[9]+d.v[5]/2-d.v[1]-d.v[8],
        d.v[10]-d.v[5]/2,
        -(d.v[5]+d.v[11]),
        d.v[3]+d.v[4]-d.v[6]-2(d.v[1]+d.v[8])),1,Val(:USQ))
end
function UnitSystem(d::Group{:USQ,<:Integer})
    Group(Values(
        -d.v[6],
        d.v[3]+d.v[4]+d.v[5]//2-d.v[1]-d.v[8],
        3(d.v[1])+2d.v[6]+4d.v[8]-d.v[3]-2d.v[4]-d.v[5]//2,
        d.v[5]//-2,
        d.v[2]+d.v[6]+d.v[7]+2(d.v[1]+d.v[8])-d.v[3]-d.v[4],
        -d.v[7],
        d.v[8],
        d.v[3]+d.v[4]+d.v[9]+d.v[5]//2-d.v[1]-d.v[8],
        d.v[10]-d.v[5]//2,
        -(d.v[5]+d.v[11]),
        d.v[3]+d.v[4]-d.v[6]-2(d.v[1]+d.v[8])),1,Val(:USQ))
end

Quantity(::ConvertUnit{D,U,S}) where {D,U,S} = Quantity{D,S}(ratio(D,U,S))
(u::UnitSystem)(::ConvertUnit{D,U,S}) where {D,U,S} = Quantity{D,u}(ratio(D,U,S))
(u::UnitSystem)(a::Number, d::Constant) = Quantity{d,u}(a)
(u::UnitSystem)(a::Number, d::AbelianGroup) = Quantity{Dimension(d),u}(a)
(u::UnitSystem)(::Constant{D}) where D = Constant{normal(u)(D)}()
(u::UnitSystem)(d::LogGroup{B}) where B = LogGroup{B}(u(d.v))
(u::UnitSystem)(d::ExpGroup{B}) where B = ExpGroup{B}(u(d.v))
