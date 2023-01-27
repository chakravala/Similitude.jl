
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

export AbelianGroup, Dimension, 𝟙, normal
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
    u∉(:permeability,:gaussgravitation) && @eval const $u = UnitSystems.$u(SI)
end

const hyperfine = SI2019(ΔνCs,inv(T))
const hubble = Hubble(𝟏,inv(T))
const cosmological = 𝟑*ΩΛ*(hubble/lightspeed(Hubble))^2
const eddington = Cosmological(𝟏,M)(QCD)
const solarmass = IAU(𝟏,M)
const earthmass = Metric(GME/G,M)(IAU)
const jupitermass = Metric(GMJ/G,M)(IAU)
const lunarmass = earthmass/μE☾
const gforce = English(𝟏,specificforce)
const atmosphere = Metric(atm,pressure)
const loschmidt = atmosphere(SI2019)/SI2019(T₀,Θ)/boltzmann(SI2019)
const amagat = loschmidt(SI2019)/avogadro(SI2019)
const wienwavelength = planck(SI)*lightspeed(SI)/boltzmann(SI)/Constant(4.965114231744276303)
const wienfrequency = Constant(2.821439372122078893)*boltzmann(SI)/planck(SI)
@pure (::typeof(loschmidt))(U::UnitSystem,P=atmosphere(U),T=SI2019(T₀,Θ)(U)) = U(P,pressure)/U(T,Θ)/boltzmann(U)
@pure mechanicalheat(U::UnitSystem) = molargas(U)*U(normal(calorie(Metric)/molargas(Metric)),Θ*N)

# angle

#const radian = Engineering(𝟏,A)
const spatian = MetricSpatian(𝟏,A)
const steradian = Engineering(𝟏,solidangle)
const degree = MetricDegree(𝟏,A)
const squaredegree = MetricDegree(𝟏,solidangle)
const gradian = MetricGradian(𝟏,A)
const bradian = Engineering(τ/𝟐^8,A)
const arcminute = MetricArcminute(𝟏,A)
const arcsecond = MetricArcsecond(𝟏,A)

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
const lunardistance = IAUE(𝟏,L)(Metric)
const jupiterdistance = IAUJ(𝟏,L)(Metric)
const parsec = astronomicalunit*(𝟐^7*𝟑^4*𝟓^3/τ)

#time

const second = Metric(𝟏,T)
const minute = (𝟐^2*𝟑*𝟓)*second
const hour = (𝟐^2*𝟑*𝟓)*minute
const day = IAU(𝟏,T)
const year = IAU(aⱼ,T)
const lightyear = year(Metric)*lightspeed(Metric)
const radarmile = 𝟐*nauticalmile(Metric)/lightspeed(Metric)
const gaussgravitation = sqrt(normal(gravitation(IAU)))*radian(IAU)/day(IAU)
const gaussianyear = turn(IAU)/gaussgravitation
const siderealyear = gaussianyear/√(solarmass+earthmass+lunarmass).v
const gaussianmonth = τ/sqrt(normal(gravitation(IAUE)))*day
const siderealmonth = gaussianmonth/normal(sqrt(earthmass(IAUE)+lunarmass(IAUE)))
const synodicmonth = inv(inv(siderealmonth(IAU))-inv(siderealyear(IAU)))
const jovianyear = τ*sqrt(normal(jupiterdistance(IAU)^3/solarmass/gravitation(IAU)))*day/normal(sqrt(solarmass+jupitermass))

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

# speed

const bubnoff = meter(Metric)/year(Metric)
const ips = IPS(𝟏,speed)
const fps = British(𝟏,speed)
const fpm = foot(British)/minute(British)
const ms = Metric(𝟏,speed)
const kmh = kilo*meter/hour
const mph = MPH(𝟏,speed)
const knot = Nautical(𝟏,speed)
const mps = mile(MPH)/second(MPH)

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
const hyl = Gravitational(𝟏,M)

# force

const dyne = Gauss(𝟏,force)
const newton = Metric(𝟏,force)
const poundal = FPS(𝟏,force)
const kilopond = Engineering(𝟏,force)
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
const talbot = Metric(𝟏,luminousenergy)
const lumerg = Gauss(centi^2*milli,luminousenergy)
const rayleigh = Metric(deka*giga,photonirradiance)
const flick = Metric(deka*giga,radiance/L)

@pure neper(U::UnitSystem) = U(𝟏,log(𝟙))
@pure bel(U::UnitSystem) = U(𝟏,log10(𝟙))
@pure decibel(U::UnitSystem) = U(𝟏,dB(𝟙))
const hertz = inv(second(Metric))
const apm = inv(minute(Metric))
const rpm = turn(Metric)/minute(Metric)
#const rpd = turn(Metric)/day(Metric)
const galileo = Gauss(𝟏,specificforce)
const eotvos = Gauss(nano,specificforce/L)
const poise = Gauss(𝟏,viscosity)
const reyn = IPS(𝟏,viscosity)
const diopter = Metric(𝟏,angularwavenumber)
const kayser = Gauss(𝟏,wavenumber)
const darcy = Gauss(milli/atm,area)
const stokes = Gauss(𝟏,diffusivity)
const katal = Metric(𝟏,catalysis)
const mpge = mile(Metric)/gasgallon(Metric)
const curie = Constant(37)*giga*hertz
const gray = Metric(𝟏,energy/M)
#const rem = centi*sievert
const roentgen = ESU(𝟏,chargedensity)(Metric)/Metric(Constant(1.293),density)
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
evaldim(::typeof(UnitSystems.solidangle)) = A^2
evaldim(unit::Group) = unit
evaldim(unit::Constant) = unit
evaldim(unit::Symbol) = evaldim(eval(unit))
evaldim(unit::Symbol,U) = evaldim(evaldim(unit),U)
evaldim(unit,U) = normal(U)(unit)

(::typeof(molarmass))(U::UnitSystem,S::UnitSystem) = (M/N)(U,S)
(::typeof(luminousefficacy))(U::UnitSystem,S::UnitSystem) = (T*J/F/L)(U,S)

dimlist(text,dim) = "$text : [$(evaldim(dim))], [$(evaldim(dim,British))], [$(evaldim(dim,Metric))], [$(evaldim(dim,EMU))], [$(evaldim(dim,ESU))]"
dimlist(U) = join(dimtext(normal(U)) == isq ? isodim.(Ref(U),usq) : unitdim.(Ref(U),usq), ", ")
isodim(U,D) = (UD = U(D); D==UD ? "$D" : "$D=$UD")
unitdim(U,D) = (io = IOBuffer(); showgroup(io,U(D),U); "$D=$(String(take!(io)))")

dimlistlatex(U) = join(isodimlatex.(Ref(U),usq), ", ")
function isodimlatex(U,D)
    UD = U(D)
    io = IOBuffer()
    latexgroup(io,D,U)
    if D≠UD
        print(io,"=")
        latexgroup(io,param(UD))
    end
    out = String(take!(io))
end

convertext(unit,fun) = """
```Julia
$(dimlist(unit,unit))
$unit(U::UnitSystem,S::UnitSystem) = $fun
$unit(v::Real,U::UnitSystem,S::UnitSystem) = v/$unit(U,S)
$(evaldim(unit)(Unified))
```
"""

@pure unitsym(x) = :nonstandard
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
$(eval(unit)(Unified))
```
"""
end

systext(sys,text) = """
```Julia
$sys = $text
$(dimlist(eval(sys)))
```
"""

# 1,2,3,4, 5, 6, 7,  8,9,10,11
#kB,ħ,𝘤,μ₀,mₑ,Mᵤ,Kcd,A,λ,αL,g₀
# F,M,L,T, Q, Θ, N,  J,A,R, C

function (u::typeof(normal(Engineering)))(d::Group)
    Group(Values(d.v[1],d.v[2],d.v[3],d.v[4],d.v[5],d.v[6],d.v[7],d.v[8],d.v[9],0,0),1,Val(:USQ))
end
function (u::typeof(normal(Gravitational)))(d::Group)
    Group(Values(d.v[1]+d.v[2],0,d.v[3]-d.v[2],d.v[4]+2(d.v[2]),d.v[5],d.v[6],d.v[7],d.v[8],0,0,0),1,Val(:USQ))
end
function (u::typeof(normal(Metric)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2],d.v[1]+d.v[3],d.v[4]-2(d.v[1]),d.v[5],d.v[6],d.v[7],d.v[8],0,0,0),1,Val(:USQ))
end
function (u::typeof(normal(MetricDegree)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2],d.v[1]+d.v[3],d.v[4]-2(d.v[1]),d.v[5],d.v[6],d.v[7],d.v[8],d.v[9],0,0),1,Val(:USQ))
end
#=function (u::typeof(normal(Astronomical)))(d::Group)
    Group(Values(d.v[1],0,d.v[3]+3d.v[2],d.v[4]-2d.v[2],d.v[5],d.v[6],d.v[7],d.v[8],d.v[9],0,0),1,Val(:USQ))
end=#

function (u::typeof(normal(Gauss)))(d::Group{<:Integer})
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]//2,d.v[1]+d.v[3]+(3//2)*d.v[5]+d.v[11],d.v[4]-2(d.v[1])-d.v[5]-d.v[11],0,d.v[6],d.v[7],d.v[8],0,0,0),1,Val(:USQ))
end
function (u::typeof(normal(ESU)))(d::Group{<:Integer})
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]//2,d.v[1]+d.v[3]+(3//2)*d.v[5],d.v[4]-2(d.v[1])-d.v[5],0,d.v[6],d.v[7],d.v[8],0,0,0),1,Val(:USQ))
end
function (u::typeof(normal(EMU)))(d::Group{<:Integer})
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]//2,d.v[1]+d.v[3]+d.v[5]//2,d.v[4]-2(d.v[1]),0,d.v[6],d.v[7],d.v[8],0,0,0),1,Val(:USQ))
end

function (u::typeof(normal(Gauss)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]//2,d.v[1]+d.v[3]+(3//2)*d.v[5]+d.v[11],d.v[4]-2(d.v[1])-d.v[5]-d.v[11],0,d.v[6],d.v[7],d.v[8],0,0,0),1,Val(:USQ))
end
function (u::typeof(normal(ESU)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]//2,d.v[1]+d.v[3]+(3//2)*d.v[5],d.v[4]-2(d.v[1])-d.v[5],0,d.v[6],d.v[7],d.v[8],0,0,0),1,Val(:USQ))
end
function (u::typeof(normal(EMU)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2]+d.v[5]//2,d.v[1]+d.v[3]+d.v[5]//2,d.v[4]-2(d.v[1]),0,d.v[6],d.v[7],d.v[8],0,0,0),1,Val(:USQ))
end

function (u::typeof(normal(Stoney)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2]+d.v[6]+d.v[7],0,d.v[3]+d.v[4]-d.v[1],d.v[5],0,0,d.v[8],0,0,0),1,Val(:USQ))
end
function (u::typeof(normal(Electronic)))(d::Group)
    Group(Values(0,0,0,d.v[3]+d.v[4]-d.v[1]-d.v[8],d.v[5],0,0,0,0,0,0),1,Val(:USQ))
end
function (u::typeof(normal(QCDoriginal)))(d::Group)
    Group(Values(0,d.v[2]+d.v[6]+d.v[7]+2(d.v[1]+d.v[8])-d.v[3]-d.v[4],0,0,d.v[5],0,0,0,0,0,0),1,Val(:USQ))
end
function (u::typeof(normal(Planck)))(d::Group)
    Group(Values(0,d.v[2]+d.v[6]+d.v[7]+2(d.v[1]+d.v[8])-d.v[3]-d.v[4],0,0,0,0,0,0,0,0,0),1,Val(:USQ))
end
function (u::typeof(normal(PlanckGauss)))(d::Group)
    Group(Values(0,d.v[2]+d.v[6]+d.v[7]+2(d.v[1]+d.v[8])-d.v[3]-d.v[4],0,0,d.v[5],0,0,0,0,0,0),1,Val(:USQ))
end
function (u::typeof(normal(Natural)))(d::Group)
    Group(Values(0,0,0,0,0,0,0,0,0,0,0),1,Val(:USQ))
end
function (u::typeof(normal(NaturalGauss)))(d::Group)
    Group(Values(0,0,0,0,d.v[5],0,0,0,0,0,0),1,Val(:USQ))
end

function (u::typeof(normal(Rydberg)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2]+d.v[7],d.v[1]+d.v[3],d.v[4]-d.v[6]+2(d.v[8]-d.v[1]),d.v[5],0,0,0,0,0,0),1,Val(:USQ))
end
function (u::typeof(normal(Hartree)))(d::Group)
    Group(Values(0,0,d.v[3]+2(d.v[4]-d.v[6])-3(d.v[1])-4(d.v[8]),0,d.v[5],0,0,0,0,0,0),1,Val(:USQ))
end
function (u::typeof(normal(Hubble)))(d::Group)
    Group(Values(0,0,0,d.v[3]+d.v[4]-d.v[1]-d.v[8],d.v[5],0,0,0,0,0,0),1,Val(:USQ))
end
function (u::typeof(normal(Cosmological)))(d::Group)
    Group(Values(0,d.v[1]+d.v[2]+d.v[6]+d.v[7],0,d.v[3]+d.v[4]-d.v[1],d.v[5],0,0,d.v[8],0,0,0),1,Val(:USQ))
end
function (u::typeof(normal(CosmologicalQuantum)))(d::Group)
    Group(Values(0,d.v[2]+d.v[6]+d.v[7]+2(d.v[1]+d.v[8])-d.v[3]-d.v[4],0,0,d.v[5],0,0,0,0,0,0),1,Val(:USQ))
end

export @unitdim, @unitgroup

"""
    @unitgroup(U::UnitSystem,S::UnitSystem) -> (u::typeof(normal(U)))(d::Group) = normal(S)(d)

Implements `Group` homomorphism for `U` in terms of existing specification from `S`.
"""
macro unitgroup(U,S)
    :((u::typeof(normal($U)))(d::Group) = normal($S)(d))
end

#=@unitgroup SI2019 Metric
@unitgroup SI1976 Metric
@unitgroup Conventional Metric
@unitgroup CODATA Metric
@unitgroup International Metric
@unitgroup InternationalMean Metric
@unitgroup MTS Metric
@unitgroup KKH Metric
@unitgroup MPH Metric
@unitgroup Nautical Metric
@unitgroup Meridian Metric
@unitgroup FFF Metric
@unitgroup IAU Metric
@unitgroup IAUE Metric
@unitgroup IAUJ Metric=#

@unitgroup MetricTurn MetricDegree
@unitgroup MetricSpatian MetricDegree
@unitgroup MetricGradian MetricDegree
@unitgroup MetricArcminute MetricDegree
@unitgroup MetricArcsecond MetricDegree
@unitgroup LorentzHeaviside Gauss
#@unitgroup Thomson EMU
#@unitgroup Kennelly EMU
@unitgroup Schrodinger Rydberg
@unitgroup QCD Planck
@unitgroup QCDGauss PlanckGauss
#@unitgroup Cosmological Hubble

#@unitgroup SI2019Engineering Engineering
#@unitgroup MeridianEngineering Engineering
#@unitgroup GravitationalSI2019 Gravitational
#@unitgroup GravitationalMeridian Gravitational
@unitgroup British Gravitational
@unitgroup English Engineering
@unitgroup Survey Engineering
@unitgroup IPS Gravitational

"""
    @unitdim(U::UnitSystem,F,M,L,T,Q,Θ,N,J="lm",A="rad")

Specify the `print` output for each base dimension of `U::UnitSystem` with `String` input arguments `force`, `mass`, `length`, `time`, `charge`, `temperature`, `molaramount`, `luminousflux`, `angle`.
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
macro unitdim(U,F,M,L,T,Q,Θ,N,J="lm",A="rad",latex=true)
    if latex
        let l = Values("\\text{$F}","\\text{$M}","\\text{$L}","\\text{$T}","\\text{$Q}","\\text{$Θ}","\\text{$N}","\\text{$J}","\\text{$A}","","")
            quote
                Similitude.dimtext(::typeof(normal($U))) = Values($F,$M,$L,$T,$Q,$Θ,$N,$J,$A,"","")
                Similitude.dimlatex(::typeof(normal($U))) = $l
            end
        end
    else
        :(Similitude.dimtext(::typeof(normal($U))) = Values($F,$M,$L,$T,$Q,$Θ,$N,$J,$A,"",""))
    end
end

@unitdim Metric "kgf" "kg" "m" "s" "C" "K" "mol"
@unitdim Meridian "kegf" "keg" "em" "s" "eC" "K" "eg-mol"
@unitdim MetricTurn "kgf" "kg" "m" "s" "C" "K" "mol" "lm" "τ" false
@unitdim MetricSpatian "kgf" "kg" "m" "s" "C" "K" "mol" "lm" "ς" false
@unitdim MetricGradian "kgf" "kg" "m" "s" "C" "K" "mol" "lm" "gon"
@unitdim MetricDegree "kgf" "kg" "m" "s" "C" "K" "mol" "lm" "deg"
@unitdim MetricArcminute "kgf" "kg" "m" "s" "C" "K" "mol" "lm" "amin"
@unitdim MetricArcsecond "kgf" "kg" "m" "s" "C" "K" "mol" "lm" "asec"
@unitdim British "lb" "slug" "ft" "s" "C" "°R" "slug-mol" "lm" "rad" false
@unitdim English "lbf" "lbm" "ft" "s" "C" "°R" "lb-mol" "lm" "rad" false
@unitdim IPS "lb" "slinch" "in" "s" "C" "°R" "slinch-mol" "lm" "rad" false
@unitdim FPS "pdl" "lb" "ft" "s" "C" "°R" "lb-mol" "lm" "rad" false
@unitdim Gauss "gf" "g" "cm" "s" "_" "K" "mol"
@unitdim IAU☉ "M☉f" "M☉" "au" "D" "C" "K" "mol" "lm" "rad" false
@unitdim IAUE "MEf" "ME" "LD" "D" "C" "K" "mol"
@unitdim IAUJ "MJf" "MJ" "JD" "D" "C" "K" "mol"
@unitdim MTS "tf" "t" "m" "s" "C" "K" "mol"
@unitdim KKH "kgf" "kg" "km" "h" "C" "K" "mol"
@unitdim MPH "lbf" "lb" "mi" "h" "C" "°R" "lb-mol" "lm" "rad" false
@unitdim Nautical "kegf" "keg" "nm" "h" "eC" "K" "eg-mol"
@unitdim FFF "firf" "fir" "fur" "ftn" "Inf" "°R" "fir-mol" "lm" "rad" false
@unitdim Hartree "F" "M" "a₀" "T" "𝘦" "Θ" "N" "J" "rad" false
@unitdim QCDoriginal "F" "mₚ" "L" "T" "𝘦" "Θ" "N" "J" "rad" false
@unitdim QCD "F" "mₚ" "L" "T" "Q" "Θ" "N" "J" "rad" false
@unitdim QCDGauss "F" "mₚ" "L" "T" "𝘦ₙ" "Θ" "N" "J" "rad" false
@unitdim PlanckGauss "F" "mP" "L" "T" "𝘦ₙ" "Θ" "N" "J" "rad" false
@unitdim NaturalGauss "F" "M" "L" "T" "𝘦ₙ" "Θ" "N" "J" "rad" false

macro unitex(U,F,M,L,T,Q,Θ,N,J="\\text{lm}",A="\\text{rad}")
    :(Similitude.dimlatex(::typeof(normal($U))) = Values($F,$M,$L,$T,$Q,$Θ,$N,$J,$A,"",""))
end

@unitex MetricTurn "\\text{kgf}" "\\text{kg}" "\\text{m}" "\\text{s}" "\\text{C}" "\\text{K}" "\\text{mol}" "\\text{lm}" "\\tau"
@unitex MetricSpatian "\\text{kgf}" "\\text{kg}" "\\text{m}" "\\text{s}" "\\text{C}" "\\text{K}" "\\text{mol}" "\\text{lm}" "\\varsigma"
@unitex British "\\text{lb}" "\\text{slug}" "\\text{ft}" "\\text{s}" "\\text{C}" "^\\circ\\text{R}" "\\text{slug-mol}"
@unitex English "\\text{lbf}" "\\text{lbm}" "\\text{ft}" "\\text{s}" "\\text{C}" "^\\circ\\text{R}" "\\text{lb-mol}"
@unitex IPS "\\text{lb}" "\\text{slinch}" "\\text{in}" "\\text{s}" "\\text{C}" "^\\circ\\text{R}" "\\text{slinch-mol}"
@unitex FPS "\\text{pdl}" "\\text{lb}" "\\text{ft}" "\\text{s}" "\\text{C}" "^\\circ\\text{R}" "\\text{lb-mol}"
@unitex IAU☉ "\\text{M}_\\odot \\text{f}" "\\text{M}_\\odot" "\\text{au}" "\\text{D}" "\\text{C}" "\\text{K}" "\\text{mol}"
@unitex MPH "\\text{lbf}" "\\text{lb}" "\\text{mi}" "\\text{h}" "\\text{C}" "^\\circ\\text{R}" "\\text{lb-mol}"
@unitex FFF "\\text{firf}" "\\text{fir}" "\\text{fur}" "\\text{ftn}" "\\infty" "^\\circ\\text{R}" "\\text{fir-mol}"
@unitex Hartree "F" "M" "\\text{a}_0" "T" "\\text{e}" "Θ" "N" "J"
@unitex QCDoriginal "F" "\\text{m}_\\text{p}" "L" "T" "\\text{e}" "Θ" "N" "J"
@unitex QCD "F" "\\text{m}_\\text{p}" "L" "T" "Q" "Θ" "N" "J"
@unitex QCDGauss "F" "\\text{m}_\\text{p}" "L" "T" "\\text{e}_\\text{n}" "Θ" "N" "J"
@unitex PlanckGauss "F" "\\text{m}_\\text{P}" "L" "T" "\\text{e}_\\text{n}" "Θ" "N" "J"
@unitex NaturalGauss "F" "M" "L" "T" "\\text{e}_\\text{n}" "Θ" "N" "J"

"""
    @unitdim(U::UnitSystem,S::UnitSystem) -> dimtext(::typeof(normal(U))) = dimtext(normal(S))

Specify the `print` output for each base dimension of `U` upon prior existing `S` data.
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
    quote
        Similitude.dimtext(::typeof(normal($U))) = dimtext(normal($S))
        Similitude.dimlatex(::typeof(normal($U))) = dimlatex(normal($S))
    end
end

@unitdim SI2019 Metric
@unitdim SI1976 Metric
@unitdim CODATA Metric
@unitdim Conventional Metric
@unitdim International Metric
@unitdim InternationalMean Metric
@unitdim Engineering Metric
@unitdim Gravitational Metric
#@unitdim SI2019Engineering Engineering
#@unitdim GravitationalSI2019 Gravitational
#@unitdim MeridianEngineering Meridian
#@unitdim GravitationalMeridian Meridian
@unitdim Survey English
@unitdim EMU Gauss
@unitdim ESU Gauss
@unitdim LorentzHeaviside Gauss
#@unitdim Thomson Gauss
#@unitdim Kennelly Metric

"""
    @unitdim(D,U,S) -> showgroup(io::IO,::typeof(U(D)),::typeof(normal(U))) = print(io,S)

Specify the `print` output `S::String` for derived dimension `D` in `U::UnitSystem`.
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
macro unitdim(D, U, S, L)
    quote
        Similitude.showgroup(io::IO,::typeof($U($D)),::typeof(normal($U))) = print(io,$S)
        Similitude.latexgroup(io::IO,::typeof($U($D)),::typeof(normal($U))) = print(io,$L)
    end
end

for U ∈ (:Engineering,:Gravitational)
    @eval begin
        @unitdim frequency $U "Hz" "\\text{Hz}"
        @unitdim frequencydrift $U "Hz⋅s⁻¹" "\\text{Hz} \\cdot \\text{s}^{-1}"
        @unitdim photonirradiance $U "Hz⋅m⁻²" "\\text{Hz} \\cdot \\text{m}^{-2}"
        @unitdim illuminance $U "lx" "\\text{lx}"
        @unitdim luminousexposure $U "lx⋅s" "\\text{lx} \\cdot \\text{s}"
        showgroup(io::IO,::typeof(luminance),::typeof(normal($U))) = print(io,"nt")
        latexgroup(io::IO,::typeof(luminance),::typeof(normal($U))) = print(io,"\\text{nt}")
    end
end
for U ∈ (:Engineering,:English,:Survey)
    @eval @unitdim specificforce $U "g₀" "\\text{g}_0"
end
for U ∈ (:Metric, :SI2019, :CODATA, :Conventional, :International, :InternationalMean, :MetricTurn, :MetricSpatian, :MetricGradian, :MetricDegree, :MetricArcminute, :MetricArcsecond)
    @eval begin
        @unitdim frequency $U "Hz" "\\text{Hz}"
        @unitdim frequencydrift $U "Hz⋅s⁻¹" "\\text{Hz} \\cdot \\text{s}^{-1}"
        @unitdim photonirradiance $U "Hz⋅m⁻²" "\\text{Hz} \\cdot \\text{m}^{-2}"
        @unitdim force $U "N" "\\text{N}"
        @unitdim inv(force) $U "N⁻¹" "\\text{N}^{-1}"
        @unitdim pressure $U "Pa" "\\text{Pa}"
        @unitdim compressibility $U "Pa⁻¹" "\\text{Pa}^{-1}"
        @unitdim energy $U "J" "\\text{J}"
        @unitdim inv(energy) $U "J⁻¹" "\\text{J}^{-1}"
        @unitdim power $U "W" "\\text{W}"
        @unitdim inv(power) $U "W⁻¹" "\\text{W}^{-1}"

        @unitdim electricpotential $U "V" "\\text{V}"
        @unitdim inv(electricpotential) $U "V⁻¹" "\\text{V}^{-1}"
        @unitdim capacitance $U "F" "\\text{F}"
        @unitdim inv(capacitance) $U "F⁻¹" "\\text{F}^{-1}"
        @unitdim resistance $U "Ω" "\\Omega"
        @unitdim conductance $U "S" "\\text{S}"
        @unitdim magneticflux $U "Wb" "\\text{Wb}"
        @unitdim inv(magneticflux) $U "Hz⋅V⁻¹" "\\text{Hz} \\cdot \\text{V}^{-1}"
        @unitdim magneticfluxdensity $U "T" "\\text{T}"
        @unitdim inv(magneticfluxdensity) $U "T⁻¹" "\\text{T}^{-1}"
        @unitdim permeance $U "H" "\\text{H}"
        @unitdim reluctance $U "H⁻¹" "\\text{H}^{-1}"

        @unitdim catalysis $U "kat" "\\text{kat}"
        @unitdim molarenergy $U "J⋅mol⁻¹" "\\text{J} \\cdot \\text{mol}^{-1}"
        @unitdim molarentropy $U "J⋅K⁻¹mol⁻¹" "\\text{J} \\cdot \\text{K}^{-1} \\text{mol}^{-1}"

        @unitdim luminousflux/power $U "lm⋅W⁻¹" "\\text{lm} \\cdot \\text{W}^{-1}"
        @unitdim power/luminousflux $U "W⋅lm⁻¹" "\\text{W} \\cdot \\text{lm}^{-1}"
        @unitdim illuminance $U "lx" "\\text{lx}"
        @unitdim luminousexposure $U "lx⋅s" "\\text{lx} \\cdot \\text{s}"

        @unitdim action*speed $U "J⋅m" "\\text{J} \\cdot \\text{m}"
        @unitdim impulse $U "N⋅s" "\\text{N} \\cdot \\text{s}"
        @unitdim yank $U "N⋅s⁻¹" "\\text{N} \\cdot \\text{s}^{-1}"
        @unitdim fluence $U "N⋅m⁻¹" "\\text{N} \\cdot \\text{m}^{-1}"
        @unitdim compliance $U "m⋅N⁻¹" "\\text{m} \\cdot \\text{N}^{-1}"

        @unitdim viscosity $U "Pa⋅s" "\\text{Pa} \\cdot \\text{s}"
        @unitdim irradiance $U "W⋅m⁻²" "\\text{W} \\cdot \\text{m}^{-2}"
        @unitdim inv(irradiance) $U "W⁻¹m²" "\\text{W}^{-1} \\text{m}^2"
        @unitdim powerdensity $U "W⋅m⁻³" "\\text{W} \\cdot \\text{m}^{-3}"
        @unitdim spectralexposure $U "J⋅m⁻²⋅Hz⁻¹" "\\text{J} \\cdot \\text{m}^{-2} \\cdot \\text{Hz}^{-1}"
        @unitdim irradiance/Θ^4 $U "W⋅m⁻²K⁻⁴" "\\text{W} \\cdot \\text{m}^{-2} \\text{K}^{-4}"
        @unitdim pressure/Θ^4 $U "J⋅m⁻³K⁻⁴" "\\text{J} \\cdot \\text{m}^{-3} \\text{K}^{-4}"
        @unitdim 𝟙/T/Θ $U "Hz⋅K⁻¹" "\\text{Hz} \\cdot \\text{K}^{-1}"
        @unitdim entropy/Q $U "V⋅K⁻¹" "\\text{V} \\cdot \\text{K}^{-1}"
        @unitdim entropy $U "J⋅K⁻¹" "\\text{J} \\cdot \\text{K}^{-1}"
        @unitdim specificentropy $U "J⋅K⁻¹kg⁻¹" "\\text{J} \\cdot \\text{K}^{-1} \\text{kg}^{-1}"
        @unitdim specificenergy $U "J⋅kg⁻¹" "\\text{J} \\cdot \\text{kg}^{-1}"
        @unitdim thermalconductivity $U "W⋅m⁻¹K⁻¹" "\\text{W} \\cdot \\text{m}^{-1} \\text{K}^{-1}"
        @unitdim thermalconductance $U "W⋅K⁻¹" "\\text{W} \\cdot \\text{K}^{-1}"
        @unitdim thermalresistance $U "K⋅W⁻¹" "\\text{K} \\cdot \\text{W}^{-1}"
        @unitdim thermalresistivity $U "K⋅m⋅W⁻¹" "\\text{K} \\cdot \\text{m} \\cdot \\text{W}^{-1}"
        @unitdim molarconductivity $U "S⋅m²mol⁻¹" "\\text{S} \\cdot \\text{m}^2 \\text{mol}^{-1}"

        @unitdim electricpotential/M $U "V⋅kg⁻¹" "\\text{V} \\cdot \\text{kg}^{-1}"
        @unitdim electricflux $U "V⋅m" "\\text{V} \\cdot \\text{m}"
        @unitdim electricfield $U "V⋅m⁻¹" "\\text{V} \\cdot \\text{m}^{-1}"
        @unitdim permittivity $U "F⋅m⁻¹" "\\text{F} \\cdot \\text{m}^{-1}"
        @unitdim inv(permittivity) $U "m⋅F⁻¹" "\\text{m} \\cdot \\text{F}^{-1}"
        @unitdim permeability $U "H⋅m⁻¹" "\\text{H} \\cdot \\text{m}^{-1}"
        @unitdim inv(permeability) $U "m⋅H⁻¹" "\\text{m} \\cdot \\text{H}^{-1}"
        @unitdim resistivity $U "Ω⋅m" "\\Omega \\cdot \\text{m}"
        @unitdim conductivity $U "S⋅m⁻¹" "\\text{S} \\cdot \\text{m}^{-1}"
        @unitdim vectorpotential $U "Wb⋅m⁻¹" "\\text{Wb} \\cdot \\text{m}^{-1}"
        @unitdim magneticmoment $U "Wb⋅m" "\\text{Wb} \\cdot \\text{m}"
        @unitdim mobility $U "m²s⁻¹V⁻¹" "\\text{m}^2 \\text{s}^{-1} \\text{V}^{-1}"
    end
end
for U ∈ (:Metric, :SI2019, :CODATA, :Conventional, :International, :InternationalMean)
    @eval begin
        @unitdim luminousintensity $U "cd" "\\text{cd}"
        Similitude.showgroup(io::IO,::typeof(luminance),::typeof(normal($U))) = print(io,"nt")
        Similitude.latexgroup(io::IO,::typeof(luminance),::typeof(normal($U))) = print(io,"\\text{nt}")
        @unitdim angularmomentum $U "J⋅s" "\\text{J} \\cdot \\text{s}"
        @unitdim magneticdipolemoment $U "J⋅T⁻¹" "\\text{J} \\cdot \\text{T}^{-1}"
    end
end
for U ∈ (:MetricTurn,:MetricSpatian,:MetricDegree,:MetricGradian,:MetricArcminute,:MetricArcsecond)
    let u = dimtext(normal(eval(U)))[9]
        let uu = dimlatex(normal(eval(U)))[9]
            @eval begin
                @unitdim angularmomentum $U $("J⋅s⋅$(u)⁻¹") $("\\text{J} \\cdot \\text{s} \\cdot $(uu)^{-1}")
                @unitdim magneticdipolemoment $U $("J⋅T⁻¹⋅$(u)⁻¹") $("\\text{J} \\cdot \\text{T}^{-1} \\cdot $(uu)^{-1}")
                @unitdim photonintensity $U $("Hz⋅$(u)⁻²") $("\\text{Hz} \\cdot $(uu)^{-2}")
                @unitdim photonradiance $U $("Hz⋅m⁻²⋅$(u)⁻²") $("\\text{Hz} \\cdot \\text{m}^{-2} \\cdot $(uu)^{-2}")
                @unitdim radiance $U $("W⋅m⁻²⋅$(u)⁻²") $("\\text{W} \\cdot \\text{m}^{-2} \\cdot $(uu)^{-2}")
                @unitdim radiance*T $U $("W⋅m⁻²⋅$(u)⁻²⋅Hz⁻¹") $("\\text{W} \\cdot \\text{m}^{-2} \\cdot $(uu)^{-2} \\cdot \\text{Hz}^{-1}")
                @unitdim radiance/L $U $("W⋅m⁻³⋅$(u)⁻²") $("\\text{W} \\cdot \\text{m}^{-3} \\cdot $(uu)^{-2}")
                @unitdim radiantintensity $U "W⋅$(u)⁻²" $("\\text{W} \\cdot $(uu)^{-2}")
                @unitdim radiantintensity*T $U $("W⋅$(u)⁻²⋅Hz⁻¹") $("\\text{W} \\cdot $(uu)^{-2} \\cdot \\text{Hz}^{-1}")
                @unitdim radiantintensity/L $U $("W⋅$(u)⁻²⋅m⁻¹") $("\\text{W} \\cdot $(uu)^{-2} \\cdot \\text{m}^{-1}")
            end
        end
    end
end

@unitdim frequency  Meridian "Hz" "\\text{Hz}"
@unitdim frequencydrift Meridian "Hz⋅s⁻¹" "\\text{Hz} \\cdot \\text{s}^{-1}"
@unitdim photonirradiance Meridian "Hz⋅m⁻²" "\\text{Hz} \\cdot \\text{m}^{-2}"
@unitdim force Meridian "eN" "\\text{eN}"
@unitdim inv(force) Meridian "eN⁻¹" "\\text{eN}^{-1}"
@unitdim pressure Meridian "ePa" "\\text{ePa}"
@unitdim compressibility Meridian "ePa⁻¹" "\\text{ePa}^{-1}"
@unitdim energy Meridian "eJ" "\\text{eJ}"
@unitdim inv(energy) Meridian "eJ⁻¹" "\\text{eJ}^{-1}"
@unitdim power Meridian "eW" "\\text{eW}"
@unitdim inv(power) Meridian "eW⁻¹" "\\text{eW}^{-1}"

@unitdim electricpotential Meridian "eV" "\\text{eV}"
@unitdim inv(electricpotential) Meridian "eV⁻¹" "\\text{eV}^{-1}"
@unitdim capacitance Meridian "eF" "\\text{eF}"
@unitdim inv(capacitance) Meridian "eF⁻¹" "\\text{eF}^{-1}"
@unitdim resistance Meridian "eΩ" "\\text{e}\\Omega"
@unitdim conductance Meridian "eS" "\\text{eS}"
@unitdim magneticflux Meridian "eWb" "\\text{eWb}"
@unitdim inv(magneticflux) Meridian "Hz⋅eV⁻¹" "\\text{Hz} \\cdot \\text{eV}^{-1}"
@unitdim magneticfluxdensity Meridian "eT" "\\text{eT}"
@unitdim inv(magneticfluxdensity) Meridian "eT⁻¹" "\\text{eT}^{-1}"
@unitdim permeance Meridian "eH" "\\text{eH}"
@unitdim reluctance Meridian "eH⁻¹" "\\text{eH}^{-1}"

@unitdim catalysis Meridian "ekat" "\\text{ekat}"
@unitdim molarenergy Meridian "eJ⋅eg-mol⁻¹" "\\text{eJ} \\cdot \\text{eg-mol}^{-1}"
@unitdim molarentropy Meridian "eJ⋅K⁻¹eg-mol⁻¹" "\\text{eJ} \\cdot \\text{K}^{-1} \\text{eg-mol}^{-1}"

@unitdim luminousflux/power Meridian "lm⋅eW⁻¹" "\\text{lm} \\cdot \\text{eW}^{-1}"
@unitdim luminousintensity Meridian "cd" "\\text{cd}"
@unitdim illuminance Meridian "elx" "\\text{elx}"
@unitdim luminousexposure Meridian "lx⋅s" "\\text{lx} \\cdot \\text{s}"
showgroup(io::IO,::typeof(luminance),::typeof(normal(Meridian))) = print(io,"ent")
latexgroup(io::IO,::typeof(luminance),::typeof(normal(Meridian))) = print(io,"\\text{ent}")

@unitdim impulse Meridian "eN⋅s" "\\text{eN} \\cdot \\text{s}"
@unitdim angularmomentum Meridian "eJ⋅s" "\\text{eJ} \\cdot \\text{s}"
@unitdim action*speed Meridian "eJ⋅em" "\\text{eJ} \\cdot \\text{em}"
@unitdim yank Meridian "eN⋅s⁻¹" "\\text{eN} \\cdot \\text{s}^{-1}"
@unitdim fluence Meridian "eN⋅em⁻¹" "\\text{eN} \\cdot \\text{em}^{-1}"
@unitdim compliance Meridian "em⋅eN⁻¹" "\\text{em} \\cdot \\text{eN}^{-1}"

@unitdim viscosity Meridian "ePa⋅s" "\\text{ePa} \\cdot \\text{s}"
@unitdim irradiance Meridian "eW⋅em⁻²" "\\text{eW} \\cdot \\text{em}^{-2}"
@unitdim inv(irradiance) Meridian "eW⁻¹em²" "\\text{eW}^{-1} \\text{em}^2"
@unitdim powerdensity Meridian "eW⋅m⁻³" "\\text{eW} \\cdot \\text{m}^{-3}"
@unitdim spectralexposure Meridian "eJ⋅em⁻²⋅Hz⁻¹" "\\text{eJ} \\cdot \\text{em}^{-2} \\cdot \\text{Hz}^{-1}"
@unitdim irradiance/Θ^4 Meridian "eW⋅em⁻²K⁻⁴" "\\text{eW} \\cdot \\text{em}^{-2} \\text{K}^{-4}"
@unitdim pressure/Θ^4 Meridian "eJ⋅em⁻³K⁻⁴" "\\text{eJ} \\cdot \\text{em}^{-3} \\text{K}^{-4}"
@unitdim 𝟙/T/Θ Meridian "Hz⋅K⁻¹" "\\text{Hz} \\cdot \\text{K}^{-1}"
@unitdim entropy/Q Meridian "eV⋅K⁻¹" "\\text{eV} \\cdot \\text{K}^{-1}"
@unitdim entropy Meridian "eJ⋅K⁻¹" "\\text{eJ} \\cdot \\text{K}^{-1}"
@unitdim specificentropy Meridian "eJ⋅K⁻¹keg⁻¹" "\\text{eJ} \\cdot \\text{K}^{-1} \\text{keg}^{-1}"
@unitdim specificenergy Meridian "eJ⋅keg⁻¹" "\\text{eJ} \\cdot \\text{keg}^{-1}"
@unitdim thermalconductivity Meridian "eW⋅em⁻¹K⁻¹" "\\text{eW} \\cdot \\text{em}^{-1} \\text{K}^{-1}"
@unitdim thermalresistance Meridian "K⋅eW⁻¹" "\\text{K} \\cdot \\text{eW}^{-1}"
@unitdim thermalresistivity Meridian "K⋅em⋅eW⁻¹" "\\text{K} \\cdot \\text{em} \\cdot \\text{eW}^{-1}"
@unitdim molarconductivity Meridian "eS⋅em²eg-mol⁻¹" "\\text{eS} \\cdot \\text{em}^2 \\text{eg-mol}^{-1}"

@unitdim electricpotential/M Meridian "eV⋅keg⁻¹" "\\text{eV} \\cdot \\text{keg}^{-1}"
@unitdim action*speed/Q Meridian "eV⋅em" "\\text{eV} \\cdot \\text{em}"
@unitdim electricfield Meridian "eV⋅em⁻¹" "\\text{eV} \\cdot \\text{em}^{-1}"
@unitdim permittivity Meridian "eF⋅em⁻¹" "\\text{eF} \\cdot \\text{em}^{-1}"
@unitdim inv(permittivity) Meridian "em⋅eF⁻¹" "\\text{em} \\cdot \\text{eF}^{-1}"
@unitdim permeability Meridian "eH⋅em⁻¹" "\\text{eH} \\cdot \\text{em}^{-1}"
@unitdim inv(permeability) Meridian "em⋅eH⁻¹" "\\text{em} \\cdot \\text{eH}^{-1}"
@unitdim resistivity Meridian "eΩ⋅em" "\\text{e}\\Omega \\cdot \\text{em}"
@unitdim conductivity Meridian "eS⋅em⁻¹" "\\text{eS} \\cdot \\text{em}^{-1}"
@unitdim magneticdipolemoment Meridian "eJ⋅eT⁻¹" "\\text{eJ} \\cdot \\text{eT}^{-1}"
@unitdim vectorpotential Meridian "eWb⋅em⁻¹" "\\text{eWb} \\cdot \\text{em}^{-1}"
@unitdim magneticmoment Meridian "eWb⋅em" "\\text{eWb} \\cdot \\text{em}"
@unitdim mobility Meridian "em²s⁻¹eV⁻¹" "\\text{em}^2 \\text{s}^{-1} \\text{eV}^{-1}"

for U ∈ (:Gauss, :EMU, :ESU, :LorentzHeaviside)
    @eval begin
        @unitdim volume $U "mL" "\\text{mL}"
        @unitdim numberdensity $U "mL⁻¹" "\\text{mL}^{-1}"
        @unitdim frequency $U "Hz" "\\text{Hz}"
        @unitdim photonirradiance $U "Hz⋅m⁻²" "\\text{Hz} \\cdot \\text{m}^{-2}"
        @unitdim force $U "dyn" "\\text{dyn}"
        @unitdim inv(force) $U "dyn⁻¹" "\\text{dyn}^{-1}"
        @unitdim specificforce $U "gal" "\\text{gal}"
        @unitdim specificforce/L $U "gal⋅cm⁻¹" "\\text{gal} \\cdot \\text{cm}^{-1}"
        @unitdim pressure $U "Ba" "\\text{Ba}"
        @unitdim compressibility $U "Ba⁻¹" "\\text{Ba}^{-1}"
        @unitdim energy $U "erg" "\\text{erg}"
        @unitdim inv(energy) $U "erg⁻¹" "\\text{erg}^{-1}"
        @unitdim power $U "erg⋅s⁻¹" "\\text{erg} \\cdot \\text{s}^{-1}"
        @unitdim inv(power) $U "s⋅erg⁻¹" "\\text{s} \\cdot \\text{erg}^{-1}"

        @unitdim catalysis $U "kat" "\\text{kat}"
        @unitdim molarenergy $U "erg⋅mol⁻¹" "\\text{erg} \\cdot \\text{mol}^{-1}"
        @unitdim molarentropy $U "erg⋅K⁻¹mol⁻¹" "\\text{erg} \\cdot \\text{K}^{-1} \\text{mol}^{-1}"

        @unitdim luminousflux/power $U "lm⋅s⋅erg⁻¹" "\\text{lm} \\cdot \\text{s} \\cdot \\text{erg}^{-1}"
        @unitdim power/luminousflux $U "erg⋅s⁻¹lm⁻¹" "\\text{erg} \\cdot \\text{s}^{-1} \\text{lm}^{-1}"
        @unitdim luminousintensity $U "cd" "\\text{cd}"
        @unitdim illuminance $U "ph" "\\text{ph}"
        showgroup(io::IO,::typeof(luminance),::typeof(normal($U))) = print(io,"sb")
        latexgroup(io::IO,::typeof(luminance),::typeof(normal($U))) = print(io,"\\text{sb}")

        @unitdim angularmomentum $U "erg⋅s" "\\text{erg} \\cdot \\text{s}"
        @unitdim action*speed $U "erg⋅cm" "\\text{erg} \\cdot \\text{cm}"
        @unitdim fluence $U "dyn⋅cm⁻¹" "\\text{dyn} \\cdot \\text{cm}^{-1}"
        @unitdim compliance $U "cm⋅dyn⁻¹" "\\text{cm} \\cdot \\text{dyn}^{-1}"
        @unitdim impulse $U "dyn⋅s" "\\text{dyn} \\cdot \\text{s}"
        @unitdim yank $U "dyn⋅s⁻¹" "\\text{dyn} \\cdot \\text{s}^{-1}"

        @unitdim viscosity $U "P" "\\text{P}"
        @unitdim diffusivity $U "St" "\\text{St}"
        @unitdim irradiance $U "erg⋅s⁻¹cm⁻²" "\\text{erg} \\cdot \\text{s}^{-1} \\text{cm}^{-2}"
        @unitdim inv(irradiance) $U "erg⁻¹s⋅cm²" "\\text{erg}^{-1} \\text{s} \\cdot \\text{cm}^2"
        @unitdim powerdensity $U "erg⋅s⁻¹mL⁻¹" "\\text{erg} \\cdot \\text{s}^{-1} \\text{mL}^{-1}"
        @unitdim spectralexposure $U "erg⋅cm⁻²⋅Hz⁻¹" "\\text{erg} \\cdot \\text{cm}^{-2} \\cdot \\text{Hz}^{-1}"
        @unitdim irradiance/Θ^4 $U "erg⋅s⁻¹cm⁻²K⁻⁴" "\\text{erg} \\cdot \\text{s}^{-1} \\text{cm}^{-2} \\text{K}^{-4}"
        @unitdim pressure/Θ^4 $U "Ba⋅K⁻⁴" "\\text{Ba} \\cdot \\text{K}^{-4}"
        @unitdim 𝟙/T/Θ $U "Hz⋅K⁻¹" "\\text{Hz} \\cdot \\text{K}^{-1}"
        @unitdim entropy $U "erg⋅K⁻¹" "\\text{erg} \\cdot \\text{K}^{-1}"
        @unitdim specificentropy $U "erg⋅K⁻¹g⁻¹" "\\text{erg} \\cdot \\text{K}^{-1} \\text{g}^{-1}"
        @unitdim specificenergy $U "erg⋅g⁻¹" "\\text{erg} \\cdot \\text{g}^{-1}"
        @unitdim thermalconductance $U "erg⋅s⁻¹K⁻¹" "\\text{erg} \\cdot \\text{s}^{-1} \\text{K}^{-1}"
        @unitdim thermalresistance $U "K⋅s⋅erg⁻¹" "\\text{K} \\cdot \\text{s} \\cdot \\text{erg}^{-1}"
        @unitdim thermalconductivity $U "erg⋅s⁻¹cm⁻¹K⁻¹" "\\text{erg} \\cdot \\text{s}^{-1} \\text{cm}^{-1} \\text{K}^{-1}"
        @unitdim thermalresistivity $U "K⋅cm⋅s⋅erg⁻¹" "\\text{K} \\cdot \\text{cm} \\cdot \\text{s} \\cdot \\text{erg}^{-1}"
    end
end

#@unitdim current EMU "Bi"
@unitdim magneticflux EMU "Mx" "\\text{Mx}"
@unitdim magneticfluxdensity EMU "G" "\\text{G}"
#@unitdim magneticfield EMU "Oe" "\\text{Oe}"
#@unitdim reluctance EMU "Bi⋅Mx⁻¹" "\\text{Bi} \\cdot \\text{Mx}^{-1}"
@unitdim magneticdipolemoment EMU "erg⋅G⁻¹" "\\text{erg} \\cdot \\text{G}^{-1}"
@unitdim vectorpotential EMU "Mx⋅cm⁻¹" "\\text{Mx} \\cdot \\text{cm}^{-1}"
#@unitdim magneticmoment EMU "Mx⋅cm" "\\text{Mx} \\cdot \\text{cm}"
#@unitdim polestrength EMU "pole" "\\text{pole}"

#@unitdim charge Gauss "Fr" "\\text{Fr}"
@unitdim magneticflux Gauss "Mx" "\\text{Mx}"
@unitdim magneticfluxdensity Gauss "G" "\\text{G}"
#@unitdim magneticfield Gauss "Oe" "\\text{Oe}"
#@unitdim reluctance Gauss "Fr⋅s⁻¹Mx⁻¹" "\\text{Fr} \\cdot \\text{s}^{-1} \\text{Mx}^{-1}"
@unitdim magneticdipolemoment Gauss "erg⋅G⁻¹" "\\text{erg} \\cdot \\text{G}^{-1}"
@unitdim vectorpotential Gauss "Mx⋅cm⁻¹" "\\text{Mx} \\cdot \\text{cm}^{-1}"
#@unitdim magneticmoment Gauss "Mx⋅cm" "\\text{Mx} \\cdot \\text{cm}"

@unitdim force MTS "sn" "\\text{sn}"
@unitdim inv(force) MTS "sn⁻¹" "\\text{sn}^{-1}"
@unitdim pressure MTS "pz" "\\text{pz}"
@unitdim compressibility MTS "pz⁻¹" "\\text{pz}^{-1}"

@unitdim mass Gravitational "hyl" "\\text{hyl}"
@unitdim mass British "slug" "\\text{slug}"
@unitdim mass IPS "slinch" "\\text{slinch}"
@unitdim molarmass Gravitational "hyl⋅mol⁻¹" "\\text{hyl} \\cdot \\text{mol}^{-1}"
@unitdim molarmass British "slug⋅slug-mol⁻¹" "\\text{slug} \\cdot \\text{slug-mol}^{-1}"
@unitdim molarmass IPS "slinch-slinch-mol⁻¹" "\\text{slinch} \\cdot \\text{slinch-mol}^{-1}"
@unitdim force FPS "pdl" "\\text{pdl}"
@unitdim pressure FPS "pdl⋅ft⁻²" "\\text{pdl} \\cdot \\text{ft}^{-2}"
@unitdim density British "slug⋅ft⁻³" "\\text{slug} \\cdot \\text{ft}^{-3}"
@unitdim density IPS "slinch⋅in⁻³" "\\text{slinch} \\cdot \\text{in}^{-3}"
@unitdim density Gravitational "hyl⋅m⁻³" "\\text{hyl} \\cdot \\text{m}^{-3}"

@unitdim L Rydberg "a₀" "\\text{a}_0"
@unitdim inv(L) Rydberg "a₀⁻¹" "\\text{a}_0^{-1}"
@unitdim area Rydberg "a₀²" "\\text{a}_0^2"
@unitdim fuelefficiency Rydberg "a₀⁻²" "\\text{a}_0^{-2}"
@unitdim volume Rydberg "a₀³" "\\text{a}_0^3"
@unitdim numberdensity Rydberg "a₀⁻³" "\\text{a}_0^{-3}"
@unitdim Q Electronic "𝘦" "\\text{e}"
@unitdim Q Stoney "𝘦" "\\text{e}"
@unitdim Q Schrodinger "𝘦" "\\text{e}"
@unitdim Q CosmologicalQuantum "𝘦ₙ" "\\text{e}_\\text{n}"
@unitdim inv(Q) Electronic "𝘦⁼¹" "\\text{e}^{-1}"
@unitdim inv(Q) Stoney "𝘦⁼¹" "\\text{e}^{-1}"
@unitdim inv(Q) Schrodinger "𝘦⁼¹" "\\text{e}^{-1}"
@unitdim inv(Q) CosmologicalQuantum "𝘦ₙ⁼¹" "\\text{e}_\\text{n}^{-1}"
@unitdim Q^2 Electronic "𝘦²" "\\text{e}^2"
@unitdim Q^2 Stoney "𝘦²" "\\text{e}^2"
@unitdim Q^2 Schrodinger "𝘦²" "\\text{e}^2"
@unitdim Q^2 CosmologicalQuantum "𝘦ₙ²" "\\text{e}_\\text{n}^2"
@unitdim inv(Q^2) Electronic "𝘦⁼²" "\\text{e}^{-2}"
@unitdim inv(Q^2) Stoney "𝘦⁼²" "\\text{e}^{-2}"
@unitdim inv(Q^2) Schrodinger "𝘦⁼²" "\\text{e}^{-2}"
@unitdim inv(Q^2) CosmologicalQuantum "𝘦ₙ⁼²" "\\text{e}_\\text{n}^{-2}"

for U ∈ (:FPS,:IPS,:British,:English,:Survey)
    @eval begin
        @unitdim frequency $U "Hz" "\\text{Hz}"
        @unitdim frequencydrift $U "Hz⋅s⁻¹" "\\text{Hz} \\cdot \\text{s}^{-1}"
        @unitdim 𝟙/T/Θ $U "Hz⋅°R⁻¹" "\\text{Hz} \\cdot ^\\circ\\text{R}^{-1}"
    end
end
for U ∈ (:FPS,:British,:English,:Survey)
    @eval begin
        @unitdim luminousintensity $U "cd" "\\text{cd}"
        @unitdim illuminance $U "fc" "\\text{fc}"
    end
end

if haskey(ENV,"UNITDOCS")
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
julia> T(MPH,Metric) # s⋅h⁻¹
$(T(MPH,Metric))

julia> T(IAU,Metric) # s⋅D⁻¹
$(T(IAU,Metric))

julia> T(Hubble,Metric)
$(T(Hubble,Metric))
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
end
