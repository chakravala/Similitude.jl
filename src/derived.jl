
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

# constants

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
#const loschmidt = atmosphere(SI2019)/SI2019(T₀,Θ)/boltzmann(SI2019)
loschmidt(U::UnitSystem,P=atmosphere(U),T=SI2019(T₀,Θ)(U)) = U(P,pressure)/U(T,Θ)/boltzmann(U)
const amagat = loschmidt(SI2019)/avogadro(SI2019)
const wienwavelength = planck(SI)*lightspeed(SI)/boltzmann(SI)/Constant(4.965114231744276303)
const wienfrequency = Constant(2.821439372122078893)*boltzmann(SI)/planck(SI)
#@pure (::typeof(loschmidt))(U::UnitSystem,P=atmosphere(U),T=SI2019(T₀,Θ)(U)) = U(P,pressure)/U(T,Θ)/boltzmann(U)
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

evaldim(::typeof(Constant(angle))) = A
evaldim(::typeof(Constant(length))) = L
evaldim(::typeof(Constant(time))) = T
evaldim(d::typeof(molarmass)) = dimensions(d)
#evaldim(::typeof(luminousefficacy)) = T*J/F/L
evaldim(::typeof(UnitSystems.solidangle)) = A^2
evaldim(::typeof(Constant(loschmidt))) = L^-3
evaldim(unit::Function) = evaldim(Constant(unit))
evaldim(unit::Group) = Constant(unit)
evaldim(unit::Constant) = unit
evaldim(unit::Symbol) = evaldim(eval(unit))
evaldim(unit::Symbol,U) = evaldim(evaldim(unit),U)
evaldim(unit::Function,U) = evaldim(evaldim(unit),U)
evaldim(unit,U) = normal(U)(unit)

(d::typeof(molarmass))(U::UnitSystem,S::UnitSystem) = dimensions(d)(U,S)
#(::typeof(luminousefficacy))(U::UnitSystem,S::UnitSystem) = (T*J/F/L)(U,S)

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
        @eval @pure unitsym(::typeof($(Constant(eval(unit))))) = $(QuoteNode(unit))
    else
        @eval @pure unitsym(::typeof($(Constant(evaldim(unit))))) = $(QuoteNode(unit))
    end
end

function unitext(unit,text)
    dim = Dimension(eval(unit))
    sym = unitsym(Constant(dim))
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

# unitdim

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
