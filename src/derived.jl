
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

export AbelianGroup, Dimension, 𝟙
export UnitSystems, Quantity, Group, LogGroup, ExpGroup
export universe, Universe, unitname, normal, dB, Db

for i ∈ 1:dims
    @eval export $(Symbol(isq[i]))
end

for u ∈ (Constants...,Physics...)
    u≠:permeability && @eval const $u = UnitSystems.$u(SI)
end

const hyperfine = SI2019(ΔνCs,inv(T))
const hubble = Hubble(𝟏,inv(T))
const cosmological = 𝟑*ΩΛ*(hubble/lightspeed(Hubble))^2
const standardgravity = Metric(g₀,L/(T*T))
const standardpressure = Metric(atm,F/(L*L))
const standardtemperature = Metric(Tₛ,Θ)
const solarmass = IAU(𝟏,M)
const earthmass = Metric(GME/G,M)(IAU)
const jupitermass = Metric(GMJ/G,M)(IAU)
const lunarmass = earthmass/μE☾
const astronomicalunit = IAU(𝟏,L)
const lunardistance = Metric(LD,L)
const mile = EE(𝟐^5*𝟑*𝟓*𝟏𝟏,L)
const clarkemile = EE(𝟐^6*𝟓*𝟏𝟗,L)
const nauticalmile = Metric(𝟐^4*𝟓^5/𝟑^3,L)
const meancalorie = InternationalMean(𝟐^2*𝟓*𝟑^2/𝟒𝟑,F*L)(Metric)
const kilocalorie = International(𝟐^5*𝟓^4*𝟑^2/𝟒𝟑,F*L)(Metric)
const calorie = kilocalorie*milli
const thermalunit = kilocalorie*ratio(M*Θ,English,Metric)
const tonsrefrigeration = thermalunit/Metric(𝟑/𝟐/𝟓,T)
const boilerhorsepower = Constant(1339)/Metric(𝟐^4*𝟑^2,T)*thermalunit
const electricalhorsepower = Metric(Constant(746),F*L/T)
const second = Metric(𝟏,T)
const minute = (𝟐^2*𝟑*𝟓)*second
const hour = (𝟐^2*𝟑*𝟓)*minute
const parsec = astronomicalunit*(𝟐^7*𝟑^4*𝟓^3/τ)
const gallon = EE(𝟕*𝟏𝟏/𝟐^6/𝟑^2,L*L*L)
const litre = Metric(milli,L*L*L)
const horsepower = British(𝟐*𝟓^2*𝟏𝟏,F*L/T)
const horsepowerwatt = British(𝟐^4*𝟑^3/𝟓*τ,F*L/T)
const horsepowermetric = GM(𝟑*𝟓^2,F*L/T)
const inchmercury = Metric(inv(inHg),F/(L*L))
const torr = Metric(atm/𝟐^3/𝟓/𝟏𝟗,F/(L*L))
const day = IAU(𝟏,T)
const year = IAU(aⱼ,T)
const gaussianyear = IAU(τ/k,T)
const siderealyear = IAU(τ/k/√(solarmass+earthmass+lunarmass).v,T)
const lightyear = year*lightspeed(IAU)

for CAL ∈ (:calₜₕ,:cal₄,:cal₁₀,:cal₂₀,:calₘ,:calᵢₜ)
    KCAL = Symbol(:k,CAL)
    @eval const $CAL = SI(UnitSystems.$CAL,energy)
    @eval const $KCAL = SI(UnitSystems.$KCAL,energy)
end

#const H0 = hubble

evaldim(::typeof(angle)) = A
evaldim(::typeof(length)) = L
evaldim(::typeof(time)) = T
evaldim(::typeof(molarmass)) = M/N
evaldim(::typeof(luminousefficacy)) = T*J/F/L
evaldim(unit::Dimension) = unit
evaldim(unit::Symbol) = evaldim(eval(unit))
evaldim(unit::Symbol,U) = evaldim(evaldim(unit),U)
evaldim(unit,U) = normal(U)(unit)

(::typeof(molarmass))(U::UnitSystem,S::UnitSystem) = (M/N)(U,S)
(::typeof(luminousefficacy))(U::UnitSystem,S::UnitSystem) = (T*J/F/L)(U,S)

convertext(unit,fun) = """
```Julia
$unit = [$(evaldim(unit))], [$(evaldim(unit,Metric))], [$(evaldim(unit,British))]
$unit(U::UnitSystem,S::UnitSystem) = $fun
$unit(v::Real,U::UnitSystem,S::UnitSystem) = v/$unit(U,S)
```
"""

# 1,2,3,4, 5, 6, 7,  8,9,10,11
#kB,ħ,𝘤,μ₀,mₑ,Mᵤ,Kcd,A,λ,αL,g₀
# F,M,L,T, Q, Θ, N,  J,A,Λ, C

function (u::typeof(normal(MetricEngineering)))(d::Group)
    Group(Values(d.v[1],d.v[2],d.v[3],d.v[4],d.v[5],d.v[6],d.v[7],d.v[8],d.v[9],0,0))
end
function (u::typeof(normal(GravitationalMetric)))(d::Group)
    Group(Values(d.v[1]+d.v[2],0,d.v[3]-d.v[2],d.v[4]+2(d.v[2]),d.v[5],d.v[6],d.v[7],d.v[8],0,0,0))
end
function (u::typeof(normal(Metric)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2],d.v[1]+d.v[3],d.v[4]-2(d.v[1]),d.v[5],d.v[6],d.v[7],d.v[8],0,0,0))
end
function (u::typeof(normal(Astronomical)))(d::Group)
    Group(Values(d.v[1],0,d.v[3]+3d.v[2],d.v[4]-2d.v[2],d.v[5],d.v[6],d.v[7],d.v[8],d.v[9],0,0))
end

function (u::typeof(normal(Gauss)))(d::Group{<:Integer})
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]//2,d.v[1]+d.v[3]+(3//2)*d.v[5]+d.v[11],d.v[4]-2(d.v[1]+d.v[5])-d.v[11],0,d.v[6],0,d.v[8],0,0,0))
end
function (u::typeof(normal(ESU)))(d::Group{<:Integer})
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]//2,d.v[1]+d.v[3]+(3//2)*d.v[5],d.v[4]-2(d.v[1]+d.v[5]),0,d.v[6],0,d.v[8],0,0,0))
end
function (u::typeof(normal(EMU)))(d::Group{<:Integer})
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]//2,d.v[1]+d.v[3]+d.v[5]//2,d.v[4]-d.v[5]-2(d.v[1]),0,d.v[6],0,d.v[8],0,0,0))
end
function (u::typeof(normal(Kennelly)))(d::Group{<:Integer})
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]//2,d.v[1]+d.v[3]+d.v[5]//2,d.v[4]-d.v[5]-2(d.v[1]),0,d.v[6],d.v[7],d.v[8],0,0,0))
end

function (u::typeof(normal(Gauss)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]/2,d.v[1]+d.v[3]+(3/2)*d.v[5]+d.v[11],d.v[4]-2(d.v[1])-d.v[5]-d.v[11],0,d.v[6],d.v[7],d.v[8],0,0,0))
end
function (u::typeof(normal(ESU)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]/2,d.v[1]+d.v[3]+(3/2)*d.v[5],d.v[4]-2(d.v[1])-d.v[5],0,d.v[6],0,d.v[8],0,0,0))
end
function (u::typeof(normal(EMU)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]/2,d.v[1]+d.v[3]+d.v[5]/2,d.v[4]-2(d.v[1]),0,d.v[6],0,d.v[8],0,0,0))
end
function (u::typeof(normal(Kennelly)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]/2,d.v[1]+d.v[3]+d.v[5]/2,d.v[4]-2(d.v[1]),0,d.v[6],d.v[7],d.v[8],0,0,0))
end

function (u::typeof(normal(Stoney)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2]+d.v[6],d.v[3]+d.v[4]-d.v[1],0,d.v[5],0,0,0,0,0,0))
end
function (u::typeof(normal(Electronic)))(d::Group)
    Group(Values(0,0,d.v[3]+d.v[4]-d.v[1],0,d.v[5],0,0,0,0,0,0))
end
function (u::typeof(normal(QCDoriginal)))(d::Group)
    Group(Values(0,d.v[2]+d.v[6]+2(d.v[1])-d.v[3]-d.v[4],0,0,d.v[5],0,0,0,d.v[9],0,0))
end
function (u::typeof(normal(PlanckGauss)))(d::Group)
    Group(Values(0,d.v[2]+d.v[6]+2(d.v[1])-d.v[3]-d.v[4],0,0,0,0,0,0,d.v[9],0,0))
end
function (u::typeof(normal(Planck)))(d::Group)
    Group(Values(0,d.v[2]+d.v[5]+d.v[6]+2(d.v[1])-d.v[3]-d.v[4],0,0,0,0,0,0,0,0,0))
end
function (u::typeof(normal(Natural)))(d::Group)
    Group(Values(0,0,0,0,0,0,0,0,0,0,0))
end
function (u::typeof(normal(NaturalGauss)))(d::Group)
    Group(Values(0,0,0,0,0,0,0,0,d.v[9],0,0))
end

function (u::typeof(normal(Rydberg)))(d::Group)
    Group(Values(0,d.v[2]+d.v[4]+d.v[6]-d.v[1],d.v[3]+2(d.v[4]-d.v[6])-3(d.v[1]),0,d.v[5],0,0,0,0,0,0))
end
function (u::typeof(normal(Hartree)))(d::Group)
    Group(Values(0,0,d.v[3]+2(d.v[4]-d.v[6])-3(d.v[1]),0,d.v[5],0,0,0,0,0,0))
end
function (u::typeof(normal(Hubble)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2],d.v[3]+d.v[4]-d.v[1],0,d.v[5],d.v[6],d.v[7],d.v[8],0,0,0))
end
function (u::typeof(normal(CosmologicalQuantum)))(d::Group)
    Group(Values(0,d.v[2]+2(d.v[1])-d.v[3]-d.v[4],0,0,d.v[5],d.v[6],d.v[7],d.v[8],0,0,0))
end

(u::typeof(normal(LorentzHeaviside)))(d::Group) = normal(Gauss)(d)
(u::typeof(normal(Thomson)))(d::Group) = normal(EMU)(d)
#(u::typeof(normal(Kennelly)))(d::Group) = normal(EMU)(d)
(u::typeof(normal(Schrodinger)))(d::Group) = normal(Rydberg)(d)
(u::typeof(normal(QCDGauss)))(d::Group) = normal(PlanckGauss)(d)
(u::typeof(normal(Cosmological)))(d::Group) = normal(Hubble)(d)

(u::typeof(normal(SI2019Engineering)))(d::Group) = normal(MetricEngineering)(d)
(u::typeof(normal(GravitationalSI2019)))(d::Group) = normal(GravitationalMetric)(d)
(u::typeof(normal(British)))(d::Group) = normal(GravitationalMetric)(d)
(u::typeof(normal(British2019)))(d::Group) = normal(GravitationalMetric)(d)
(u::typeof(normal(English)))(d::Group) = normal(MetricEngineering)(d)
(u::typeof(normal(English2019)))(d::Group) = normal(MetricEngineering)(d)
(u::typeof(normal(Survey)))(d::Group) = normal(MetricEngineering)(d)
(u::typeof(normal(Survey2019)))(d::Group) = normal(MetricEngineering)(d)

@doc """
$(convertext(:length,"planck(U,S)/mass(U,S)/speed(U,S)"))

Extent of one-dimensional shape or `length` (m), unit conversion factor.

```Julia
julia> L(CGS,Metric) # m⋅cm⁻¹
$(L(CGS,Metric))

julia> L(IAU,Metric) # m⋅au⁻¹
$(L(IAU,Metric))

julia> L(English,Metric) # m⋅ft⁻¹
$(L(English,Metric))

julia> L(EnglishUS,English) # ft⋅ftUS⁻¹
$(L(EnglishUS,English))

julia> L(PlanckGauss,Metric) # m⋅ℓP⁻¹
$(L(PlanckGauss,Metric))
```
""" L, ft, ftUS

@doc """
$(convertext(:time,"length(U,S)/speed(U,S)"))

Dimension along which events are ordered or `T` (s), unit conversion factor.

```Julia
julia> T(IAU,Metric) # s⋅day⁻¹
$(T(IAU,Metric))

julia> T(PlanckGauss,Metric) # s⋅tP⁻¹
$(T(PlanckGauss,Metric))
```
""" T
