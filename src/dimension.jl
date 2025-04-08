
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
@pure dimension(D) = D

value(::Constant{D}) where D = value(D)

showgroup(io::IO,x::Constant{D},u::UnitSystem) where D = showgroup(io,D,dimtext(u),'𝟙')
showgroup(io::IO,x::AbelianGroup,u::UnitSystem) = showgroup(io,Constant(x),u)
showgroup(io::IO,x::Group,u::UnitSystem) = showgroup(io,Constant(x),u)
showgroup(io::IO,x,u,c) = showgroup(io,x,u)

import FieldAlgebra: latexgroup
latexgroup(io::IO,x::Constant{D},u) where D = latexgroup(io,D,dimlatex(u),"\\mathbb{1}")
latexgroup(io::IO,x::AbelianGroup,u::UnitSystem) = latexgroup(io,Constant(x),u)
latexgroup(io::IO,x::Group,u::UnitSystem) = latexgroup(io,Constant(x),u)
latexgroup(io::IO,x,u,c) = latexgroup(io,x,u)

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
Base.:+(a::Constant,b::Group{:Constants}) = param(a)+b
Base.:+(a::Group{:Constants},b::Constant) = a+param(b)
Base.:-(a::Constant,b::Group{:Constants}) = param(a)-b
Base.:-(a::Group{:Constants},b::Constant) = a-param(b)

Base.:*(a::Group{:Constants},b::Group{:USQ}) = Group(b.v,a*b.c,Val(:USQ))
Base.:*(a::Group{:USQ},b::Group{:Constants}) = Group(a.v,a.c*b,Val(:USQ))
Base.:/(a::Group{:Constants},b::Group{:USQ}) = a*inv(b)
Base.:/(a::Group{:USQ},b::Group{:Constants}) = a*inv(b)

@group USQ F M L T Q Θ N J A R C

const isq = Values('F','M','L','T','Q','Θ','N','J','A','R','C')
const dims = length(isq)

const dimensionless = if CONSTDIM
    Constant(valueat(0,dims,:USQ))
else
    valueat(0,dims,:USQ)
end
const 𝟙 = dimensionless
const USQ = Values(F,M,L,T,Q,Θ,N,J,A,R,C)
const usq = USQ
const usqlatex = Values("\\text{k}_\\text{B}","\\text{N}_\\text{A}","\\hbar","\\text{c}","\\text{e}","\\text{K}_\\text{cd}","\\Delta\\nu_\\text{Cs}","\\text{R}_{\\infty}","\\alpha","\\mu_\\text{eu}","\\mu_\\text{pu}","\\Omega_{\\Lambda}","\\text{H}_0","\\text{g}_0","\\text{a}_\\text{j}","\\text{au}","\\text{ft}","\\text{ft}_\\text{US}","\\text{lb}","\\text{T}_0","\\text{atm}","\\text{in}_\\text{Hg}","\\text{R}_\\text{K}^{90}","\\text{K}_\\text{J}^{90}","\\text{R}_\\text{K}","\\text{K}_\\text{J}","\\text{R}_\\text{u}","\\Omega_\\text{it}","\\text{V}_\\text{it}","\\text{k}_\\text{G}","\\text{m}_\\text{P}","\\text{GM}_\\text{E}","\\text{GM}_\\text{J}","\\varphi","\\gamma","e","\\tau","2","3","5","7","11","19","43")

naturalunits(U) = [x=>x(U) for x ∈ Similitude.USQ]
morphism(U) = [getproperty.(param.(U.(Similitude.usq)),:v)[j][i] for i ∈ 1:11, j ∈ 1:11]

#FieldAlgebra.latext(::Group{:USQ}) = Values('F','M','L','Q',"\\Theta",'N','J','A','R','C')
FieldAlgebra.latext(::Group{:Constants}) = usqlatex
@group2 Constants begin
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

if CONSTVAL
    phys(j,k=vals) = Constant(valueat(j,k,:Constants))
else
    phys(j,k=vals) = valueat(j,k,:Constants)
end

const basis = Values("kB", "NA", "𝘩", "𝘤", "𝘦", "Kcd", "ΔνCs", "R∞", "α", "μₑᵤ", "μₚᵤ", "ΩΛ", "H0", "g₀", "aⱼ", "au", "ft", "ftUS", "lb", "T₀", "atm", "inHg", "RK90", "KJ90", "RK", "KJ", "Rᵤ2014", "Ωᵢₜ", "Vᵢₜ", "kG", "mP", "GME", "GMJ", "φ", "γ", "ℯ", "τ", "2", "3", "5", "7", "11", "19","43")
const vals = length(basis)

# convertunit

struct ConvertUnit{U,S,D} <: AbstractModule
    v::D
    @pure ConvertUnit{U,S}(v::D) where {U,S,D} = new{normal(U),normal(S),D}(v)
end

dimensions(c::ConvertUnit) = c.v
convertdim(c::ConvertUnit{U,S}) where {U,S} = convertdim(dimensions(c),U,S)
convertdim(::Constant{D},U,S) where D = convertdim(D,U,S)
dimconvert(x::T,d,u,s) where T = (isone(ratio(d,u,s)) ? zero(x) : x)::T
convertdim(d::Group{:USQ,T},U,S) where T = Constant{Group{:USQ,T}(dimconvert.(d.v,usq,Ref(U),Ref(S)))}()

function Base.show(io::IO,c::ConvertUnit{U,S}) where {U,S}
    d = convertdim(dimensions(c),U,S)
    print(io, ratio(dimensions(c),U,S), " [")
    showgroup(io,S(d),S)
    print(io, "]/[")
    showgroup(io,U(d),U)
    print(io, "] ", unitname(U), " -> ", unitname(S))
end

@pure Base.inv(::ConvertUnit{U,S}) where {U,S} = ConvertUnit{U,S}(inv(dimensions(c)))

(D::Constant)(U::UnitSystem,S::UnitSystem) = ConvertUnit{U,S}(D)
(D::Constant)(U::UnitSystem) = U(ratio(D,Natural,U),D)
(D::Constant)(v::Real,U::UnitSystem,S::UnitSystem=Metric) = v/ratio(D,U,S) #U(_,D)
(D::Group{:USQ})(U::UnitSystem,S::UnitSystem) = ConvertUnit{U,S}(D)
(D::Group{:USQ})(U::UnitSystem) = U(ratio(D,Natural,U),D)
(D::Group{:USQ})(v::Real,U::UnitSystem,S::UnitSystem=Metric) = v/ratio(D,U,S) #U(_,D)
#(D::Dimension)(v::Number,U::UnitSystem,S::UnitSystem=Natural) = U(v/ratio(D,U,S),D)

Base.log(x::ConvertUnit{U,S}) where {U,S} = ConvertUnit{U,S}(log(dimensions(x)))
Base.log2(x::ConvertUnit{U,S}) where {U,S} = ConvertUnit{U,S}(log2(dimensions(x)))
Base.log10(x::ConvertUnit{U,S}) where {U,S} = ConvertUnit{U,S}(log10(dimensions(x)))
Base.log(b::Number,x::ConvertUnit{U,S}) where {U,S} = ConvertUnit{U,S}(log(b,dimensions(x)))
Base.exp(x::ConvertUnit{U,S}) where {U,S} = ConvertUnit{U,S}(exp(dimensions(x)))
Base.exp2(x::ConvertUnit{U,S}) where {U,S} = ConvertUnit{U,S}(exp2(dimensions(x)))
Base.exp10(x::ConvertUnit{U,S}) where {U,S} = ConvertUnit{U,S}(exp10(dimensions(x)))
Base.:^(a::Number,b::ConvertUnit{U,S}) where {U,S} = ConvertUnit{U,S}(a^dimensions(b))
Base.:^(a::ConvertUnit{U,S},b::Integer) where {U,S} = Quantity{U,S}(dimensions(a)^b)
Base.:^(a::ConvertUnit{U,S},b::Rational{Int}) where {U,S} = Quantity{U,S}(dimensions(a)^b)
Base.:*(a::ConvertUnit{U,S},b::ConvertUnit{U,S}) where {U,S} = ConvertUnit{U,S}(dimensions(a)*dimensions(b))
Base.:/(a::ConvertUnit{U,S},b::ConvertUnit{U,S}) where {U,S} = ConvertUnit{U,S}(dimensions(a)/dimensions(b))

# quantity

struct Quantity{U,T,D} <: AbstractModule
    v::T
    d::D
    Quantity{U,T}(v,d::D) where {T,U,D} = new{normal(U),T,D}(v,d)
    @pure Quantity{U,Int}(v,d::D) where {U,D} = new{normal(U),Int,D}(v,d)
    @pure Quantity{U,Float64}(v,d::D) where {U,D} = new{normal(U),Float64,D}(v,d)
end

Quantity{U}(v::Quantity,d) where U = Quantity{U}(v.v,d)
Quantity{U}(v::T,d) where {U,T} = Quantity{U,T}(v,d)
@pure Quantity{U}(v::Int,d::D) where {U,D} = Quantity{U,Int}(v,d)::Quantity{normal(U),Int,D}
@pure Quantity{U}(v::Float64,d::D) where {U,D} = Quantity{U,Float64}(v,d)::Quantity{normal(U),Float64,D}
@pure Quantity(U::UnitSystem,v::Int,d) = Quantity{U}(v,d)
@pure Quantity(U::UnitSystem,v::Float64,d) = Quantity{U}(v,d)
Quantity(U::UnitSystem,v,d) = Quantity{U}(v,d)
@pure Quantity(v::Int,d,U::UnitSystem) = Quantity{U}(v,d)
@pure Quantity(v::Float64,d,U::UnitSystem) = Quantity{U}(v,d)
Quantity(v,d,U::UnitSystem) = Quantity{U}(v,d)
@pure Quantity(d,U::UnitSystem,v::Int) = Quantity{U}(v,d) # deprecate
@pure Quantity(d,U::UnitSystem,v::Float64) = Quantity{U}(v,d) # deprecate
Quantity(d,U::UnitSystem,v) = Quantity{U}(v,d) # deprecate

Base.convert(::Type{Float64},x::Quantity) = convert(Float64,x.v)
Base.convert(::Type{Float64},x::Quantity{U,Float64}) where U = x.v
Base.convert(::Type{T},x::Quantity{U,T}) where {U,T<:Number} = x.v

import UnitSystems: evaldim
UnitSystems.isquantity(U::Quantity) = true

quantity(q) = q
quantity(q::Quantity) = q.v
@pure quantity(q::Quantity{U,Int} where U) = q.v
@pure quantity(q::Quantity{U,Float64} where U) = q.v
@pure dimensions(q::Quantity) = q.d
@pure Dimension(q::Quantity) = q.d
@pure Dimension(q) = q
@pure unitsystem(::Quantity{U}) where U = dimension(U)
@pure unitsystem2(::Quantity{U}) where U = U

function Base.show(io::IO,x::Quantity{U}) where U
    print(io, x.v, " [")
    showgroup(io,normal(U)(dimensions(x)),U)
    print(io, "] ", unitname(U))
end

Base.log(x::Quantity{U}) where U = Quantity{U}(log(x.v),log(dimensions(x)))
Base.log2(x::Quantity{U}) where U = Quantity{U}(log2(x.v),log2(dimensions(x)))
Base.log10(x::Quantity{U}) where U = Quantity{U}(log10(x.v),log10(dimensions(x)))
Base.log(b::Number,x::Quantity{U}) where U = Quantity{U}(log(b,x.v),log(b,dimensions(x)))
Base.exp(x::Quantity{U}) where U = Quantity{U}(exp(x.v),exp(dimensions(x)))
Base.exp2(x::Quantity{U}) where U = Quantity{U}(exp2(x.v),exp2(dimensions(x)))
Base.exp10(x::Quantity{U}) where U = Quantity{U}(exp10(x.v),exp10(dimensions(x)))
#Base.:^(a::Constant,b::Quantity{U}) where U = Quantity{U}(a^b.v,a^dimensions(x))
Base.:^(a::UnitSystems.Constant,b::Quantity{U}) where U = Quantity{U}(a^b.v,a^dimensions(b))
Base.:^(a::Number,b::Quantity{U}) where U = Quantity{U}(a^b.v,a^dimensions(b))
Base.:^(a::Quantity{U},b::Number) where U = Quantity{U}(a.v^b,dimensions(a)^b)
Base.:^(a::Quantity{U},b::Integer) where U = Quantity{U}(a.v^b,dimensions(a)^b)
Base.:^(a::Quantity{U},b::Rational{Int}) where U = Quantity{U}(a.v^b,dimensions(a)^b)
Base.:+(a::Number,b::Quantity{U}) where U = U(D)==𝟙 ? Quantity{U}(a+b.v) : throw(error("$(U(D)) ≠ 𝟙 "),dimensions(b))
Base.:+(a::Quantity{U},b::Number) where U = U(D)==𝟙 ? Quantity{U}(a.v+b) : throw(error("$(U(D)) ≠ 𝟙 "),dimensions(a))
Base.:-(a::Number,b::Quantity{U}) where U = U(D)==𝟙 ? Quantity{U}(a-b.v) : throw(error("$(U(D)) ≠ 𝟙 "),dimensions(a))
Base.:-(a::Quantity{U},b::Number) where U = U(D)==𝟙 ? Quantity{U}(a.v-b) : throw(error("$(U(D)) ≠ 𝟙 "),dimensions(a))
Base.:*(a::Real,b::Quantity{U}) where U = Quantity{U}(a*b.v,dimensions(b))
Base.:*(a::Quantity{U},b::Real) where U = Quantity{U}(a.v*b,dimensions(a))
Base.:*(a::Complex,b::Quantity{U}) where U = Quantity{U}(a*b.v,dimensions(b))
Base.:*(a::Quantity{U},b::Complex) where U = Quantity{U}(a.v*b,dimensions(a))
Base.:*(a::Quantity{U},b::Quantity{U}) where U = Quantity{U}(a.v*b.v,dimensions(a)*dimensions(b))
#Base.:/(a::Quantity{U},b::Quantity{U}) where U = Quantity{U}(a.v/b.v,dimensions(a)/dimensions(b))
Base.:/(a::Quantity{U},b::Quantity{U}) where U = Quantity{U}(a.v/b.v,dimensions(a)/dimensions(b))
Base.:/(a::Quantity{A},b::Quantity{B}) where {A,B} = ConvertUnit{A,B}((a.v/b.v)*dimensions(a))
Base.:/(a::Number,b::Quantity) = a*inv(b)
Base.:/(a::Quantity,b::Number) = a*inv(b)
Base.:/(a::Quantity,b::UnitSystems.Constant) = a*inv(b)
Base.:/(a::UnitSystems.Constant,b::Quantity) = a*inv(b)
#Base.:/(a::Quantity{U},b::ConvertUnit{S,U}) where {U,S} = a*ConvertUnit{U,S}(dimensions(a))
Base.:-(a::Quantity{U}) where U = Quantity{U}(-a.v,dimensions(a))
Base.inv(a::Quantity{U}) where U = Quantity{U}(inv(a.v),inv(dimensions(a)))
Base.sqrt(a::Quantity{U}) where U = Quantity{U}(sqrt(a.v),sqrt(dimensions(a)))
Base.cbrt(a::Quantity{U}) where U = Quantity{U}(cbrt(a.v),cbrt(dimensions(a)))

#Base.:*(a::Dimension{A},b::Quantity{U}) where {U,A} = Quantity{a,U}(10^-(isunknown(A) ? A.v.v[end] : A.v[end]))*b
#Base.:*(a::Quantity{U},b::Dimension{B}) where {U,B} = a*Quantity{b,U}(10^-(isunknown(B) ? B.v.v[end] : B.v[end]))
Base.:*(a::Quantity{U},b::ConvertUnit{U,U}) where U = dimensions(a)==dimensions(b) || dimensions(b)==𝟙 ? a : error("ConvertUnit incompatible $(dimensions(a)) ≠ $(dimensions(b))")
Base.:*(a::ConvertUnit{U,U},b::Quantity{U}) where U = dimensions(b)==dimensions(a) || dimensions(a)==𝟙 ? b : error("ConvertUnit incompatible $(dimensions(b)) ≠ $(dimensions(a))")
function Base.:*(a::ConvertUnit{S,U},b::Quantity{U}) where {U,S}
    A,D = dimensions(a),dimensions(b)
    (A==D && S==U) ? b : inv(A)==D ? ConvertUnit{U,S}(inv(A))*b : error("ConvertUnit incompatible $(inv(A)) ≠ $D")
end
function Base.:*(a::Quantity{U},b::ConvertUnit{S,U}) where {U,S}
    D,A = dimensions(a),dimensions(b)
    (D==A && S==U) ? a : inv(A)==D ? a*ConvertUnit{U,S}(inv(A)) : error("ConvertUnit incompatible $D ≠ $(inv(A))")
end
function Base.:*(a::ConvertUnit{U,S},b::Quantity{U}) where {U,S}
    A,D = dimensions(a),dimensions(b)
    A==D ? Quantity{S}(ratio(D,U,S)*b.v,D) : error("ConvertUnit incompatible $A ≠ $B")
end
function Base.:*(a::Quantity{U},b::ConvertUnit{U,S}) where {U,S}
    D,A = dimensions(a),dimensions(b)
    A==D ? Quantity{S}(a.v*ratio(D,U,S),D) : error("ConvertUnit incompatible $A ≠ $B")
end
function Base.:+(a::Constant,b::Quantity{U}) where U
    D = dimensions(b)
    U(D)==𝟙 ? Quantity{U}(a+b.v,D) : throw(error("$(U(D)) ≠ 𝟙 "))
end
function Base.:+(a::Quantity{U},b::Constant) where U
    D = dimensions(a)
    U(D)==𝟙 ? Quantity{U}(a.v+b,D) : throw(error("$(U(D)) ≠ 𝟙 "))
end
function Base.:-(a::Constant,b::Quantity{U}) where U
    D = dimensions(b)
    U(D)==𝟙 ? Quantity{D,U}(a-b.v) : throw(error("$(U(D)) ≠ 𝟙 "))
end
function Base.:-(a::Quantity{U},b::Constant) where U
    D = dimensions(a)
    U(D)==𝟙 ? Quantity{U}(a.v-b,D) : throw(error("$(U(D)) ≠ 𝟙 "))
end
Base.:*(a::Constant,b::Quantity{U}) where U = Quantity{U}(a*b.v,dimensions(b))
Base.:*(a::Quantity{U},b::Constant) where U = Quantity{U}(a.v*b,dimensions(a))
Base.:*(a::Group,b::Quantity{U}) where U = Quantity{U}(a*b.v,dimensions(b))
Base.:*(a::Quantity{U},b::Group) where U = Quantity{U}(a.v*b,dimensions(a))
#Base.:/(a::Constant,b::Quantity{U}) where U = Quantity{U}(a/b.v,inv(dimensions(b)))
#Base.:/(a::Quantity{U},b::Constant) where U = Quantity{U}(a.v/b,dimensions(a))

Base.:(==)(a::Quantity{U},b::Quantity{U}) where U = U(dimensions(a)) == U(dimensions(b)) && a.v == b.v
function Base.:(==)(a::Number,b::Quantity)
    D = dimensions(b)
    iszero(norm(first(value(D),dims-1))) && a == b.v*10^last(value(D))
end
function Base.:(==)(a::Quantity,b::Number)
    D = dimensions(a)
    iszero(norm(first(value(D),dims-1))) && b == a.v*10^last(value(D))
end

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
function Base.:+(a::Quantity{U},b::Quantity{U}) where U
    A,B = dimensions(a),dimensions(b)
    Quantity{U}(a.v+b.v,A≠B ? add(A,B,U) : A)
end
function Base.:-(a::Quantity{U},b::Quantity{U}) where U
    A,B = dimensions(a),dimensions(b)
    Quantity{U}(a.v-b.v,A≠B ? sub(A,B,U) : A)
end

# other

Base.isone(x::Quantity{D}) where D = false

# Quantities

struct Quantities{U,N,T,D} <: TupleVector{N,T}
    v::Values{N,T}
    d::D
    Quantities{U,N,T}(v::Values{N,T},d::D) where {D,U,N,T} = new{normal(U),N,T,D}(v,d)
end

Quantities{U,N,T}(v::Tuple,d::D) where {D,U,N,T} = Quantities{U,N,T}(Values{N,T}(v),d)
Quantities{U}(v::Values{N,T},d::D) where {D,U,N,T} = Quantities{U,N,T}(v,d)
Quantities(d,U,v) = Quantities{U}(v,d)

function (S::UnitSystem)(x::Quantities{U}) where U
    D = dimensions(x)
    Quantities{S}(D.(x.v,Ref(U),Ref(S)),D)
end

export Quantity, Quantities

Base.getindex(x::Quantities{U},i::Integer) where U = Quantity{U}(x.v[i],dimensions(x))
Base.getindex(x::Quantities{U},i::Int) where U = Quantity{U}(x.v[i],dimensions(x))

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

function Quantity(c::ConvertUnit{U,S}) where {U,S}
    D = dimensions(c)
    Quantity{S}(ratio(D,U,S),D)
end
function (u::UnitSystem)(c::ConvertUnit{U,S}) where {U,S}
    D = dimensions(c)
    Quantity{u}(ratio(D,U,S),D)
end
(u::UnitSystem)(a::Number, d::Constant) = Quantity{u}(a,d)
(u::UnitSystem)(a::Number, d::AbelianGroup) = Quantity{u}(a,d)
(u::UnitSystem)(::Constant{D}) where D = Constant{normal(u)(D)}()
#(u::UnitSystem)(d::Group) = normal(Metric)(d)
(u::UnitSystem)(d::Group) = normal(u)≠u ? normal(u)(d) : normal(Metric)(d)
(u::UnitSystem)(d::LogGroup{B}) where B = LogGroup{B}(u(d.v))
(u::UnitSystem)(d::ExpGroup{B}) where B = ExpGroup{B}(u(d.v))
