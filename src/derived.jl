
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
export universe, Universe, unitname, normal, logdb, expdb, dB

for i ∈ 1:dims
    @eval export $(Symbol(isq[i]))
end

for unit ∈ Dimensionless
    @eval @pure $unit(C::Coupling) = UnitSystems.$unit(C)
    @eval @pure $unit(U::UnitSystem) = UnitSystems.$unit(universe(U))
end

for u ∈ (Constants...,Physics...)
    u≠:permeability && @eval const $u = UnitSystems.$u(SI)
end

const hyperfine = SI2019(ΔνCs,inv(T))
const hubble = Hubble(𝟏,inv(T))
const cosmological = 𝟑*ΩΛ*(hubble/lightspeed(Hubble))^2
const solarmass = IAU(𝟏,M)
const earthmass = Metric(GME/G,M)(IAU)
const jupitermass = Metric(GMJ/G,M)(IAU)
const lunarmass = earthmass/μE☾
const gaussianyear = IAU(τ/k,T)
const siderealyear = IAU(τ/k/√(solarmass+earthmass+lunarmass).v,T)
const gforce = English(𝟏,specificforce)
const atmosphere = Metric(atm,pressure)
const loschmidt = atmosphere(SI2019)/SI2019(T₀,Θ)/boltzmann(SI2019)
const amagat = loschmidt(SI2019)/avogadro(SI2019)
const wienwavelength = planck(SI)*lightspeed(SI)/boltzmann(SI)/Constant(4.965114231744276303)
const wienfrequency = Constant(2.821439372122078893)*boltzmann(SI)/planck(SI)
@pure (::typeof(loschmidt))(U::UnitSystem,P=atmosphere(U),T=SI2019(T₀,Θ)(U)) = U(P,pressure)/U(T,Θ)/boltzmann(U)
@pure mechanicalheat(U::UnitSystem) = molargas(U)*U(normal(calorie(Metric)/molargas(Metric)),Θ*N)

# angle

const radian = MetricEngineering(𝟏,A)
const steradian = MetricEngineering(𝟏,solidangle)
const degree = MetricEngineering(τ/𝟐^3/𝟑^2/𝟓,A)
const gradian = MetricEngineering(τ/𝟐^4/𝟓^2,A)
const arcminute = degree/𝟐^2/𝟑/𝟓
const arcsecond = arcminute/𝟐^2/𝟑/𝟓

# length

const meter = Metric(𝟏,L)
const angstrom = hecto*pico*meter
const inch = IPS(𝟏,L)
#const rackunit = foot*𝟕/𝟐^4/𝟑
const foot = English(𝟏,L)
const surveyfoot = Survey(𝟏,L)
const yard = 𝟑*foot
const mile = English(𝟐^5*𝟑*𝟓*𝟏𝟏,L)
const statutemile = Survey(𝟐^5*𝟑*𝟓*𝟏𝟏,L)
const earthradius = sqrt(earthmass(Metric)*gravitation(Metric)/gforce(Metric))
const greatcircle = τ*earthradius
const earthmeter = Meridian(𝟏,L)
const nauticalmile = Nautical(𝟏,L)
const admiraltymile = English(𝟐^6*𝟓*𝟏𝟗,L)
const meridianmile = Metric(𝟐^4*𝟓^5/𝟑^3,L)
const astronomicalunit = IAU(𝟏,L)
const lunardistance = Metric(LD,L)
const parsec = astronomicalunit*(𝟐^7*𝟑^4*𝟓^3/τ)

#time

const second = Metric(𝟏,T)
const minute = (𝟐^2*𝟑*𝟓)*second
const hour = (𝟐^2*𝟑*𝟓)*minute
const day = IAU(𝟏,T)
const year = IAU(aⱼ,T)
const lightyear = year*lightspeed(IAU)
const radarmile = 𝟐*nauticalmile(Metric)/lightspeed(Metric)

# area

const barn = Metric((𝟐*𝟓)^-28,area)
const hectare = Metric(hecto*hecto,area)
const acre = MPH(𝟐^-7/𝟓,area)(English)
const surveyacre = Survey(𝟐^3*𝟑^2*𝟓*𝟏𝟏^2,area)
#const township = 𝟐^9*𝟑^2*𝟓*surveyacre
#const footballfield = English(𝟐^8*𝟑^2*𝟓^2,area)

# volume

const gallon = IPS(𝟑*𝟕*𝟏𝟏,volume)
const liter = Metric(milli,volume)
const quart = gallon/𝟐^2
const pint = quart/𝟐
const cup = pint/𝟐
const fluidounce = cup/𝟐^3
const teaspoon = 𝟓*milli*liter
const tablespoon = 𝟑*teaspoon
#const oilbarrel = 𝟐*𝟑*𝟕*gallon

# mass

const gram = Metric(milli,M)
const earthgram = Meridian(milli,M)
const kilogram = Metric(𝟏,M)
const tonne = Metric(kilo,M)
const ton = English(𝟐*kilo,M)
const pound = English(𝟏,M)
const ounce = English(𝟐^-4,M)
const grain = milli*pound/𝟕
const slug = British(𝟏,M)
const slinch = IPS(𝟏,M)
const hyl = GravitationalMetric(𝟏,M)

# force

const dyne = Gauss(𝟏,force)
const newton = Metric(𝟏,force)
const poundal = FPS(𝟏,force)
const kilopond = MetricEngineering(𝟏,force)
const poundforce = English(𝟏,F)

# pressure

const psi = IPS(𝟏,pressure)
const bar = Metric(hecto*kilo,pressure)
const barye = Gauss(𝟏,pressure)
const pascal = Metric(𝟏,pressure)
const technicalatmosphere = kilopond/(centi*meter(ME))^2
const inchmercury = Metric(inv(inHg),pressure)
const torr = Metric(atm/𝟐^3/𝟓/𝟏𝟗,pressure)

# energy

const erg = Gauss(𝟏,energy)
const joule = Metric(𝟏,energy)
const footpound = poundforce*foot
const meancalorie = InternationalMean(𝟐^2*𝟓*𝟑^2/𝟒𝟑,energy)(Metric)
const kilocalorie = International(𝟐^5*𝟓^4*𝟑^2/𝟒𝟑,energy)(Metric)
const calorie = kilocalorie*milli
const earthcalorie = mechanicalheat(Meridian)
const thermalunit = mechanicalheat(English)
const tontnt = giga*calorie(Metric)
const gasgallon = 𝟐*𝟑*𝟏𝟗*kilo*thermalunit(Metric)

# power

const watt = Metric(𝟏,power)
const horsepower = British(𝟐*𝟓^2*𝟏𝟏,power)
const horsepowerwatt = British(𝟐^4*𝟑^3/𝟓*τ,power)
const horsepowermetric = GM(𝟑*𝟓^2,power)
const tonsrefrigeration = thermalunit(Metric)/Metric(𝟑/𝟐/𝟓,T)
const boilerhorsepower = Constant(1339)/Metric(𝟐^4*𝟑^2,T)*thermalunit(Metric)
const electricalhorsepower = Metric(Constant(746),power)

# electromagnetic

const coulomb = Metric(𝟏,Q)
const ampere = Metric(𝟏,current)
const volt = Metric(𝟏,electricpotential)
const henry = Metric(𝟏,inductance)
const ohm = Metric(𝟏,resistance)
const siemens = Metric(𝟏,conductance)
const farad = Metric(𝟏,capacitance)
const weber = Metric(𝟏,magneticflux)
const tesla = Metric(𝟏,magneticfluxdensity)
const statcoulomb = ESU(𝟏,Q)
const statampere = ESU(𝟏,current)
const statvolt = ESU(𝟏,electricpotential)
const stathenry = ESU(𝟏,inductance)
const statohm = ESU(𝟏,resistance)
const statmho = ESU(𝟏,conductance)
const statfarad = ESU(𝟏,capacitance)
const statweber = ESU(𝟏,magneticflux)
const stattesla = ESU(𝟏,magneticfluxdensity)
const abcoulomb = EMU(𝟏,Q)
const abampere = EMU(𝟏,current)
const abvolt = EMU(𝟏,electricpotential)
const abhenry = EMU(𝟏,inductance)
const abohm = EMU(𝟏,resistance)
const abmho = EMU(𝟏,conductance)
const abfarad = EMU(𝟏,capacitance)
const maxwell = EMU(𝟏,magneticflux)
const gauss = EMU(𝟏,magneticfluxdensity)
const oersted = EMU(𝟏,magneticfield)
const gilbert = EMU(𝟏/𝟐/τ,current/A)
const earthcoulomb = Meridian(𝟏,Q)
const electronvolt = elementarycharge(SI2019)*SI2019(𝟏,electricpotential)

# temperature

#const freezing = Metric(T₀-milli,Θ)
const boiling = Metric(T₀+Constant(99.9839),Θ)
const sealevel = Metric(T₀+𝟑*𝟓,Θ)
const kelvin = Metric(𝟏,Θ)
const celsius = Metric(T₀,Θ)
const rankine = English(𝟏,Θ)
const fahrenheit = English(Constant(459.67),Θ)
#const delisle = Metric(𝟐/𝟑,Θ)
#const reaumur = Metric(𝟓/𝟐^2,Θ)

# mole

const mole = Metric(𝟏,N)
const earthmole = Meridian(𝟏,N)
const poundmole = English(𝟏,N)
const slugmole = British(𝟏,N)
const slinchmole = IPS(𝟏,N)

# photometric

const lumen = Metric(𝟏,luminousflux)
const candela = Metric(𝟏,luminousintensity)
const lux = Metric(𝟏,illuminance)
const phot = Gauss(𝟏,illuminance)
const footcandle = English(𝟏,illuminance)
const nit = Metric(𝟏,luminance)
const apostilb = Metric(𝟐/τ,luminance)
const stilb = Gauss(𝟏,luminance)
const lambert = Gauss(𝟐/τ,luminance)
const footlambert = English(𝟐/τ,luminance)
const bril = centi*nano*lambert

@pure neper(U::UnitSystem) = U(𝟏,log(𝟙))
@pure bel(U::UnitSystem) = U(𝟏,log10(𝟙))
@pure decibel(U::UnitSystem) = U(𝟏,dB(𝟙))
const hertz = inv(second(Metric))
const rpm = inv(minute(Metric))
#const rpd = turn(Metric)/day(Metric)
const galileo = Gauss(𝟏,specificforce)
const eotvos = Gauss(nano,specificforce/L)
const poise = Gauss(𝟏,viscosity)
const reyn = IPS(𝟏,viscosity)
const diopter = Metric(𝟏,wavenumber)
const kayser = Gauss(𝟏,wavenumber)
const darcy = Gauss(milli/atm,area)
const stokes = Gauss(𝟏,diffusivity)
const katal = Metric(𝟏,catalysis)
const mpge = mile(Metric)/gasgallon(Metric)
const curie = Constant(37)*giga*hertz
const sievert = Metric(𝟏,energy/M)
#const rem = centi*sievert
const roentgen = ESU(𝟏,chargedensity)(Metric)/Metric(Constant(1.293),density)
const bubnoff = meter(Metric)/year(Metric)
const rayl = Metric(𝟏,specificimpedance)
const langley = calorie(Metric)/(centi*meter(Metric))^2
const jansky = Metric((𝟐*𝟓)^-26,fluence)
const solarflux = hecto*hecto*jansky

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

dimlist(text,dim) = "$text : [$(evaldim(dim))], [$(evaldim(dim,British))], [$(evaldim(dim,Metric))], [$(evaldim(dim,EMU))], [$(evaldim(dim,ESU))]"

convertext(unit,fun) = """
```Julia
$(dimlist(unit,unit))
$unit(U::UnitSystem,S::UnitSystem) = $fun
$unit(v::Real,U::UnitSystem,S::UnitSystem) = v/$unit(U,S)
```
"""

@pure unitsym(x) = :nonstandard
@pure unitsym(::typeof(A)) = :angle
for unit ∈ Convert
    if unit ∉ (:length,:time,:angle,:molarmass,:luminousefficacy)
        @eval @pure unitsym(::typeof($(eval(unit)))) = $(QuoteNode(unit))
    else
        @eval @pure unitsym(::typeof($(evaldim(unit)))) = $(QuoteNode(unit))
    end
end

function unitext(unit,text)
    dim = Dimension(eval(unit))
    sym = unitsym(dim)
    return """
```Julia
$unit(U::UnitSystem) = $text
$(dimlist(sym,dim))
```
"""
end

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
#=function (u::typeof(normal(Astronomical)))(d::Group)
    Group(Values(d.v[1],0,d.v[3]+3d.v[2],d.v[4]-2d.v[2],d.v[5],d.v[6],d.v[7],d.v[8],d.v[9],0,0))
end=#

function (u::typeof(normal(Gauss)))(d::Group{<:Integer})
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]//2,d.v[1]+d.v[3]+(3//2)*d.v[5]+d.v[11],d.v[4]-2(d.v[1]+d.v[5])-d.v[11],0,d.v[6],d.v[7],d.v[8],0,0,0))
end
function (u::typeof(normal(ESU)))(d::Group{<:Integer})
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]//2,d.v[1]+d.v[3]+(3//2)*d.v[5],d.v[4]-2(d.v[1]+d.v[5]),0,d.v[6],d.v[7],d.v[8],0,0,0))
end
function (u::typeof(normal(EMU)))(d::Group{<:Integer})
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]//2,d.v[1]+d.v[3]+d.v[5]//2,d.v[4]-d.v[5]-2(d.v[1]),0,d.v[6],d.v[7],d.v[8],0,0,0))
end

function (u::typeof(normal(Gauss)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]/2,d.v[1]+d.v[3]+(3/2)*d.v[5]+d.v[11],d.v[4]-2(d.v[1])-d.v[5]-d.v[11],0,d.v[6],d.v[7],d.v[8],0,0,0))
end
function (u::typeof(normal(ESU)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]/2,d.v[1]+d.v[3]+(3/2)*d.v[5],d.v[4]-2(d.v[1])-d.v[5],0,d.v[6],d.v[7],d.v[8],0,0,0))
end
function (u::typeof(normal(EMU)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]/2,d.v[1]+d.v[3]+d.v[5]/2,d.v[4]-2(d.v[1]),0,d.v[6],d.v[7],d.v[8],0,0,0))
end

function (u::typeof(normal(Stoney)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2]+d.v[6],0,d.v[3]+d.v[4]-d.v[1],d.v[5],0,0,0,0,0,0))
end
function (u::typeof(normal(Electronic)))(d::Group)
    Group(Values(0,0,0,d.v[3]+d.v[4]-d.v[1],d.v[5],0,0,0,0,0,0))
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
    Group(Values(0,d.v[2]+d.v[4]+d.v[6]-d.v[1],0,d.v[3]+2(d.v[4]-d.v[6])-3(d.v[1]),d.v[5],0,0,0,0,0,0))
end
function (u::typeof(normal(Hartree)))(d::Group)
    Group(Values(0,0,0,d.v[3]+2(d.v[4]-d.v[6])-3(d.v[1]),d.v[5],0,0,0,0,0,0))
end
function (u::typeof(normal(Hubble)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2],0,d.v[3]+d.v[4]-d.v[1],d.v[5],d.v[6],d.v[7],d.v[8],0,0,0))
end
function (u::typeof(normal(CosmologicalQuantum)))(d::Group)
    Group(Values(0,d.v[2]+2(d.v[1])-d.v[3]-d.v[4],0,0,d.v[5],d.v[6],d.v[7],d.v[8],0,0,0))
end

export @unitdim, @unitgroup

"""
    @unitgroup(U::UnitSystem,S::UnitSystem) -> (u::typeof(normal(U)))(d::Group) = normal(S)(d)

Implements `Group` homomorphism for `U` in terms of existing specification from `S`.
"""
macro unitgroup(U,S)
    :((u::typeof(normal($U)))(d::Group) = normal($S)(d))
end

@unitgroup LorentzHeaviside Gauss
#@unitgroup Thomson EMU
@unitgroup Kennelly EMU
@unitgroup Schrodinger Rydberg
@unitgroup QCD Planck
@unitgroup QCDGauss PlanckGauss
@unitgroup Cosmological Hubble

@unitgroup SI2019Engineering MetricEngineering
@unitgroup MeridianEngineering MetricEngineering
@unitgroup GravitationalSI2019 GravitationalMetric
@unitgroup GravitationalMeridian GravitationalMetric
@unitgroup British GravitationalMetric
@unitgroup English MetricEngineering
@unitgroup Survey MetricEngineering
@unitgroup IPS GravitationalMetric

"""
    @unitdim(U::UnitSystem,F,M,L,T,Q,Θ,N,J="lm",A="rad")

Specify the `print` output for each base `Dimension` of `U::UnitSystem` with `String` input arguments `force`, `mass`, `length`, `time`, `charge`, `temperature`, `molaramount`, `luminousflux`, `angle`.
```Julia
@unitdim Gauss "gf" "g" "cm" "s" "C" "K" "mol"
@unitdim Metric "kgf" "kg" "m" "s" "C" "K" "mol"
@unitdim British "lb" "slug" "ft" "s" "C" "°R" "slug-mol"
@unitdim IPS "lb" "slinch" "in" "s" "C" "°R" "slinch-mol"
@unitdim FPS "pdl" "lb" "ft" "s" "C" "°R" "lb-mol"
@unitdim English "lbf" "lbm" "ft" "s" "C" "°R" "lb-mol"
@unitdim IAU☉ "M☉f" "M☉" "au" "D" "C" "K" "mol"
```
These standard examples are some of the built-in defaults.
"""
macro unitdim(U,F,M,L,T,Q,Θ,N,J="lm",A="rad",Λ="",C="")
    :(dimtext(::typeof(normal($U))) = Values($F,$M,$L,$T,$Q,$Θ,$N,$J,$A,$Λ,$C))
end

@unitdim Metric "kgf" "kg" "m" "s" "C" "K" "mol"
@unitdim Meridian "kegf" "keg" "em" "s" "eC" "K" "eg-mol"
@unitdim British "lb" "slug" "ft" "s" "C" "°R" "slug-mol"
@unitdim English "lbf" "lbm" "ft" "s" "C" "°R" "lb-mol"
@unitdim IPS "lb" "slinch" "in" "s" "C" "°R" "slinch-mol"
@unitdim FPS "pdl" "lb" "ft" "s" "C" "°R" "lb-mol"
@unitdim Gauss "gf" "g" "cm" "s" "_" "K" "mol"
@unitdim IAU☉ "M☉f" "M☉" "au" "D" "C" "K" "mol"
@unitdim IAUE "MEf" "ME" "au" "D" "C" "K" "mol"
@unitdim IAUJ "MJf" "MJ" "au" "D" "C" "K" "mol"
@unitdim MTS "tf" "t" "m" "s" "C" "K" "mol"
@unitdim KKH "kgf" "kg" "km" "h" "C" "K" "mol"
@unitdim MPH "lbf" "lb" "mi" "h" "C" "°R" "lb-mol"
@unitdim Nautical "kegf" "keg" "nm" "h" "eC" "K" "eg-mol"
@unitdim FFF "firf" "fir" "fur" "ftn" "Inf" "°R" "fir-mol"

"""
    @unitdim(U::UnitSystem,S::UnitSystem) -> dimtext(::typeof(normal(U))) = dimtext(normal(S))

Specify the `print` output for each base `Dimension` of `U` upon prior existing `S` data.
```Julia
@unitdim EMU Gauss
@unitdim ESU Gauss
@unitdim LorentzHeaviside Gauss
@unitdim SI2019 Metric
@unitdim SI1976 Metric
@unitdim CODATA Metric
@unitdim Conventional Metric
@unitdim International Metric
@unitdim InternationalMean Metric
@unitdim Survey English
```
These standard examples are some of the built-in defaults.
"""
macro unitdim(U,S)
    :(dimtext(::typeof(normal($U))) = dimtext(normal($S)))
end

@unitdim SI2019 Metric
@unitdim SI1976 Metric
@unitdim CODATA Metric
@unitdim Conventional Metric
@unitdim International Metric
@unitdim InternationalMean Metric
@unitdim MetricEngineering Metric
@unitdim GravitationalMetric Metric
@unitdim SI2019Engineering MetricEngineering
@unitdim GravitationalSI2019 GravitationalMetric
@unitdim MeridianEngineering Meridian
@unitdim GravitationalMeridian Meridian
@unitdim Survey English
@unitdim EMU Gauss
@unitdim ESU Gauss
@unitdim LorentzHeaviside Gauss
#@unitdim Thomson Gauss
@unitdim Kennelly Metric

"""
    @unitdim(D,U,S) -> showgroup(io::IO,::typeof(U(D)),::typeof(normal(U))) = print(io,S)

Specify the `print` output `S::String` for derived `D::Dimension` in `U::UnitSystem`.
```Julia
@unitdim magneticflux Gauss "Mx"
@unitdim magneticfluxdensity Gauss "G"
@unitdim magneticfield Gauss "Oe"
@unitdim frequency Metric "Hz"
@unitdim force Metric "N"
@unitdim pressure Metric "Pa"
@unitdim energy Metric "J"
@unitdim power Metric "W"
@unitdim mass British "slug"
@unitdim force FPS "pdl"
```
These standard examples are some of the built-in defaults.
"""
macro unitdim(D, U, S)
    :(showgroup(io::IO,::typeof($U($D)),::typeof(normal($U))) = print(io,$S))
end

for U ∈ (:MetricEngineering, :SI2019Engineering,:GravitationalMetric,:GravitationalSI2019)
    @eval begin
        @unitdim frequency $U "Hz"
        @unitdim frequencydrift $U "Hz*s⁻¹"
        @unitdim illuminance $U "lx"
        @unitdim luminousexposure $U "lx*s"
        showgroup(io::IO,::typeof(luminance),::typeof(normal($U))) = print(io,"nt")
    end
end
for U ∈ (:MetricEngineering,:SI2019Engineering,:MeridianEngineering,:English,:Survey)
    @eval @unitdim specificforce $U "g₀"
end
for U ∈ (:Metric, :SI2019, :CODATA, :Conventional, :International, :InternationalMean)
    @eval begin
        @unitdim frequency $U "Hz"
        @unitdim frequencydrift $U "Hz*s⁻¹"
        @unitdim force $U "N"
        @unitdim inv(force) $U "N⁻¹"
        @unitdim pressure $U "Pa"
        @unitdim compressibility $U "Pa⁻¹"
        @unitdim energy $U "J"
        @unitdim inv(energy) $U "J⁻¹"
        @unitdim power $U "W"
        @unitdim inv(power) $U "W⁻¹"

        @unitdim electricpotential $U "V"
        @unitdim inv(electricpotential) $U "V⁻¹"
        @unitdim capacitance $U "F"
        @unitdim inv(capacitance) $U "F⁻¹"
        @unitdim resistance $U "Ω"
        @unitdim conductance $U "S"
        @unitdim magneticflux $U "Wb"
        @unitdim inv(magneticflux) $U "Hz*V⁻¹"
        @unitdim magneticfluxdensity $U "T"
        @unitdim inv(magneticfluxdensity) $U "T⁻¹"
        @unitdim permeance $U "H"
        @unitdim reluctance $U "H⁻¹"

        @unitdim catalysis $U "kat"
        @unitdim molarenergy $U "J*mol⁻¹"
        @unitdim molarentropy $U "J*K⁻¹mol⁻¹"

        @unitdim luminousflux/power $U "lm*W⁻¹"
        @unitdim power/luminousflux $U "W*lm⁻¹"
        @unitdim luminousintensity $U "cd"
        @unitdim illuminance $U "lx"
        showgroup(io::IO,::typeof(luminance),::typeof(normal($U))) = print(io,"nt")
        @unitdim luminousexposure $U "lx*s"

        @unitdim angularmomentum $U "J*s"
        @unitdim action*speed $U "J*m"
        @unitdim impulse $U "N*s"
        @unitdim yank $U "N*s⁻¹"
        @unitdim fluence $U "N*m⁻¹"
        @unitdim compliance $U "m*N⁻¹"

        @unitdim viscosity $U "Pa*s"
        @unitdim intensity $U "W*m⁻²"
        @unitdim powerdensity $U "W*m⁻³"
        @unitdim intensity/Θ^4 $U "W*m⁻²K⁻⁴"
        @unitdim pressure/Θ^4 $U "J*m⁻³K⁻⁴"
        @unitdim 𝟙/T/Θ $U "Hz*K⁻¹"
        @unitdim entropy/Q $U "V*K⁻¹"
        @unitdim entropy $U "J*K⁻¹"
        @unitdim specificentropy $U "J*K⁻¹kg⁻¹"
        @unitdim specificenergy $U "J*kg⁻¹"
        @unitdim thermalconductivity $U "W*m⁻¹K⁻¹"
        @unitdim thermalconductance $U "W*K⁻¹"
        @unitdim thermalresistance $U "K*W⁻¹"
        @unitdim thermalresistivity $U "K*m*W⁻¹"
        @unitdim molarconductivity $U "S*m²mol⁻¹"

        @unitdim electricpotential/M $U "V*kg⁻¹"
        @unitdim electricpotential*L $U "V*m"
        @unitdim electricfield $U "V*m⁻¹"
        @unitdim permittivity $U "F*m⁻¹"
        @unitdim inv(permittivity) $U "m*F⁻¹"
        @unitdim permeability $U "H*m⁻¹"
        @unitdim inv(permeability) $U "m*H⁻¹"
        @unitdim resistivity $U "Ω*m"
        @unitdim conductivity $U "S*m⁻¹"
        @unitdim magneticdipolemoment $U "J*T⁻¹"
        @unitdim vectorpotential $U "Wb*m⁻¹"
        @unitdim magneticmoment $U "Wb*m"
        @unitdim mobility $U "m²s⁻¹V⁻¹"
    end
end

@unitdim frequency  Meridian "Hz"
@unitdim frequencydrift Meridian "Hz*s⁻¹"
@unitdim force Meridian "eN"
@unitdim inv(force) Meridian "eN⁻¹"
@unitdim pressure Meridian "ePa"
@unitdim compressibility Meridian "ePa⁻¹"
@unitdim energy Meridian "eJ"
@unitdim inv(energy) Meridian "eJ⁻¹"
@unitdim power Meridian "eW"
@unitdim inv(power) Meridian "eW⁻¹"

@unitdim electricpotential Meridian "eV"
@unitdim inv(electricpotential) Meridian "eV⁻¹"
@unitdim capacitance Meridian "eF"
@unitdim inv(capacitance) Meridian "eF⁻¹"
@unitdim resistance Meridian "eΩ"
@unitdim conductance Meridian "eS"
@unitdim magneticflux Meridian "eWb"
@unitdim inv(magneticflux) Meridian "Hz*eV⁻¹"
@unitdim magneticfluxdensity Meridian "eT"
@unitdim inv(magneticfluxdensity) Meridian "eT⁻¹"
@unitdim permeance Meridian "eH"
@unitdim reluctance Meridian "eH⁻¹"

@unitdim catalysis Meridian "ekat"
@unitdim molarenergy Meridian "eJ*eg-mol⁻¹"
@unitdim molarentropy Meridian "eJ*K⁻¹eg-mol⁻¹"

@unitdim luminousflux/power Meridian "lm*eW⁻¹"
@unitdim luminousintensity Meridian "cd"
@unitdim illuminance Meridian "elx"
@unitdim luminousexposure Meridian "lx*s"
showgroup(io::IO,::typeof(luminance),::typeof(normal(Meridian))) = print(io,"ent")

@unitdim impulse Meridian "eN*s"
@unitdim angularmomentum Meridian "eJ*s"
@unitdim action*speed Meridian "eJ*em"
@unitdim yank Meridian "eN*s⁻¹"
@unitdim fluence Meridian "eN*em⁻¹"
@unitdim compliance Meridian "em*eN⁻¹"

@unitdim viscosity Meridian "ePa*s"
@unitdim intensity Meridian "eW*em⁻²"
@unitdim powerdensity Meridian "eW*m⁻³"
@unitdim intensity/Θ^4 Meridian "eW*em⁻²K⁻⁴"
@unitdim pressure/Θ^4 Meridian "eJ*em⁻³K⁻⁴"
@unitdim 𝟙/T/Θ Meridian "Hz*K⁻¹"
@unitdim entropy/Q Meridian "eV*K⁻¹"
@unitdim entropy Meridian "eJ*K⁻¹"
@unitdim specificentropy Meridian "eJ*K⁻¹keg⁻¹"
@unitdim specificenergy Meridian "eJ*keg⁻¹"
@unitdim thermalconductivity Meridian "eW*em⁻¹K⁻¹"
@unitdim thermalresistance Meridian "K*eW⁻¹"
@unitdim thermalresistivity Meridian "K*em*eW⁻¹"
@unitdim molarconductivity Meridian "eS*em²eg-mol⁻¹"

@unitdim electricpotential/M Meridian "eV*kg⁻¹"
@unitdim action*speed/Q Meridian "eV*em"
@unitdim electricfield Meridian "eV*em⁻¹"
@unitdim permittivity Meridian "eF*em⁻¹"
@unitdim inv(permittivity) Meridian "em*eF⁻¹"
@unitdim permeability Meridian "eH*em⁻¹"
@unitdim inv(permeability) Meridian "em*eH⁻¹"
@unitdim resistivity Meridian "eΩ*em"
@unitdim conductivity Meridian "eS*em⁻¹"
@unitdim magneticdipolemoment Meridian "eJ*eT⁻¹"
@unitdim vectorpotential Meridian "eWb*em⁻¹"
@unitdim magneticmoment Meridian "eWb*em"
@unitdim mobility Meridian "em²s⁻¹eV⁻¹"

for U ∈ (:Gauss, :EMU, :ESU, :LorentzHeaviside)
    @eval begin
        @unitdim frequency $U "Hz"
        @unitdim force $U "dyn"
        @unitdim inv(force) $U "dyn⁻¹"
        @unitdim specificforce $U "gal"
        @unitdim specificforce/L $U "gal*cm⁻¹"
        @unitdim pressure $U "Ba"
        @unitdim compressibility $U "Ba⁻¹"
        @unitdim energy $U "erg"
        @unitdim inv(energy) $U "erg⁻¹"
        @unitdim power $U "erg*s⁻¹"
        @unitdim inv(power) $U "s*erg⁻¹"

        @unitdim catalysis $U "kat"
        @unitdim molarenergy $U "erg*mol⁻¹"
        @unitdim molarentropy $U "erg*K⁻¹mol⁻¹"

        @unitdim luminousflux/power $U "lm*s*erg⁻¹"
        @unitdim power/luminousflux $U "erg*s⁻¹lm⁻¹"
        @unitdim luminousintensity $U "cd"
        @unitdim illuminance $U "ph"
        showgroup(io::IO,::typeof(luminance),::typeof(normal($U))) = print(io,"sb")

        @unitdim angularmomentum $U "erg*s"
        @unitdim action*speed $U "erg*cm"
        @unitdim fluence $U "dyn*cm⁻¹"
        @unitdim compliance $U "cm*dyn⁻¹"
        @unitdim impulse $U "dyn*s"
        @unitdim yank $U "dyn*s⁻¹"

        @unitdim viscosity $U "P"
        @unitdim diffusivity $U "St"
        @unitdim intensity $U "erg*s⁻¹cm⁻²"
        @unitdim powerdensity $U "erg*s⁻¹cm⁻³"
        @unitdim intensity/Θ^4 $U "erg*s⁻¹cm⁻²K⁻⁴"
        @unitdim pressure/Θ^4 $U "Ba*K⁻⁴"
        @unitdim 𝟙/T/Θ $U "Hz*K⁻¹"
        @unitdim entropy $U "erg*K⁻¹"
        @unitdim specificentropy $U "erg*K⁻¹g⁻¹"
        @unitdim specificenergy $U "erg*g⁻¹"
        @unitdim thermalconductance $U "erg*s⁻¹K⁻¹"
        @unitdim thermalresistance $U "K*s*erg⁻¹"
        @unitdim thermalconductivity $U "erg*s⁻¹cm⁻¹K⁻¹"
        @unitdim thermalresistivity $U "K*cm*s*erg⁻¹"
    end
end

@unitdim current EMU "Bi"
@unitdim magneticflux EMU "Mx"
@unitdim magneticfluxdensity EMU "G"
@unitdim magneticfield EMU "Oe"
@unitdim reluctance EMU "Bi*Mx⁻¹"
@unitdim magneticdipolemoment EMU "erg*G⁻¹"
@unitdim vectorpotential EMU "Mx*cm⁻¹"
@unitdim magneticmoment EMU "Mx*cm"
@unitdim polestrength EMU "pole"

@unitdim charge Gauss "Fr"
@unitdim magneticflux Gauss "Mx"
@unitdim magneticfluxdensity Gauss "G"
@unitdim magneticfield Gauss "Oe"
@unitdim reluctance Gauss "Fr*s⁻¹Mx⁻¹"
@unitdim magneticdipolemoment Gauss "erg*G⁻¹"
@unitdim vectorpotential Gauss "Mx*cm⁻¹"
@unitdim magneticmoment Gauss "Mx*cm"

@unitdim force MTS "sn"
@unitdim inv(force) MTS "sn⁻¹"
@unitdim pressure MTS "pz"
@unitdim compressibility MTS "pz⁻¹"

@unitdim mass GravitationalMetric "hyl"
@unitdim mass GravitationalSI2019 "hyl"
@unitdim mass GravitationalMeridian "ehyl"
@unitdim mass British "slug"
@unitdim mass IPS "slinch"
@unitdim force FPS "pdl"
@unitdim pressure FPS "pdl*ft⁻²"
@unitdim density British "slug*ft⁻³"
@unitdim density IPS "slinch*in⁻³"
@unitdim density GravitationalMetric "hyl*m⁻³"
@unitdim density GravitationalSI2019 "hyl*m⁻³"
@unitdim density GravitationalMeridian "ehyl*m⁻³"

for U ∈ (:FPS,:IPS,:British,:English,:Survey)
    @eval begin
        @unitdim frequency $U "Hz"
        @unitdim frequencydrift $U "Hz*s⁻¹"
        @unitdim 𝟙/T/Θ $U "Hz*°R⁻¹"
    end
end
for U ∈ (:FPS,:British,:English,:Survey)
    @eval begin
        @unitdim luminousintensity $U "cd"
        @unitdim illuminance $U "fc"
    end
end

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

@doc """
    neper(U::UnitSystem) = U(𝟏,log(𝟙))

Logarithmic unit expressing the ratio of a dimensional quanty.
```Julia
julia> neper(Metric)
$(neper(Metric))

julia> exp(neper(Metric))
$(exp(neper(Metric)))
```
""" neper

@doc """
    bel(U::UnitSystem) = U(𝟏,log10(𝟙))

Logarithmic unit expressing the ratio of a dimensional quanty.
```Julia
julia> bel(Metric)
$(bel(Metric))

julia> exp10(bel(Metric))
$(exp10(bel(Metric)))
```
""" bel

@doc """
    decibel(U::UnitSystem) = U(𝟏,logdb(𝟙))

Logarithmic unit expressing the ratio of a dimensional quanty.
```Julia
julia> decibel(Metric)
$(decibel(Metric))

julia> expdb(decibel(Metric))
$(expdb(decibel(Metric)))
```
""" decibel
