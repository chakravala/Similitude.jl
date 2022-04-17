
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
const foot = English(𝟏,L)
const inch = IPS(𝟏,L)
#const rackunit = foot*𝟕/𝟐^4/𝟑
const yard = 𝟑*foot
const surveyfoot = Survey(𝟏,L)
const statutemile = Survey(𝟐^5*𝟑*𝟓*𝟏𝟏,L)
const earthradius = sqrt(earthmass(Metric)*gravitation(Metric)/gforce(Metric))
const greatcircle = τ*earthradius
const earthmeter = Meridian(𝟏,L)
const nauticalmile = Nautical(𝟏,L)
const astronomicalunit = IAU(𝟏,L)
const lunardistance = Metric(LD,L)
const mile = English(𝟐^5*𝟑*𝟓*𝟏𝟏,L)
const admiraltymile = English(𝟐^6*𝟓*𝟏𝟗,L)
const meridianmile = Metric(𝟐^4*𝟓^5/𝟑^3,L)
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
const acre = MPH(𝟐^-7/𝟓,L^2)(English)
const surveyacre = Survey(𝟐^3*𝟑^2*𝟓*𝟏𝟏^2,L^2)
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

const sealevel = Metric(T₀+𝟑*𝟓,Θ)
const kelvin = Metric(𝟏,Θ)
const rankine = English(𝟏,Θ)
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
const lambert = CGS(𝟐/τ,luminance)
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

(u::typeof(normal(LorentzHeaviside)))(d::Group) = normal(Gauss)(d)
#(u::typeof(normal(Thomson)))(d::Group) = normal(EMU)(d)
#(u::typeof(normal(Kennelly)))(d::Group) = normal(EMU)(d)
(u::typeof(normal(Schrodinger)))(d::Group) = normal(Rydberg)(d)
(u::typeof(normal(QCD)))(d::Group) = normal(Planck)(d)
(u::typeof(normal(QCDGauss)))(d::Group) = normal(PlanckGauss)(d)
(u::typeof(normal(Cosmological)))(d::Group) = normal(Hubble)(d)

(u::typeof(normal(SI2019Engineering)))(d::Group) = normal(MetricEngineering)(d)
(u::typeof(normal(MeridianEngineering)))(d::Group) = normal(MetricEngineering)(d)
(u::typeof(normal(GravitationalSI2019)))(d::Group) = normal(GravitationalMetric)(d)
(u::typeof(normal(GravitationalMeridian)))(d::Group) = normal(GravitationalMetric)(d)
(u::typeof(normal(British)))(d::Group) = normal(GravitationalMetric)(d)
#(u::typeof(normal(British2019)))(d::Group) = normal(GravitationalMetric)(d)
(u::typeof(normal(English)))(d::Group) = normal(MetricEngineering)(d)
#(u::typeof(normal(English2019)))(d::Group) = normal(MetricEngineering)(d)
(u::typeof(normal(Survey)))(d::Group) = normal(MetricEngineering)(d)
#(u::typeof(normal(Survey2019)))(d::Group) = normal(MetricEngineering)(d)
(u::typeof(normal(IPS)))(d::Group) = normal(GravitationalMetric)(d)
#(u::typeof(normal(IPS2019)))(d::Group) = normal(GravitationalMetric)(d)

dimtext(::typeof(normal(Metric))) = Values("kgf","kg","m","s","C","K","mol","lm","rad","","")
dimtext(::typeof(normal(MetricEngineering))) = dimtext(normal(Metric))
dimtext(::typeof(normal(GravitationalMetric))) = dimtext(normal(Metric))
dimtext(::typeof(normal(Meridian))) = Values("kegf","keg","em","s","eC","K","eg-mol","lm","rad","","")
dimtext(::typeof(normal(British))) = Values("lb","slug","ft","s","C","°R","slug-mol","lm","rad","","")
dimtext(::typeof(normal(English))) = Values("lbf","lbm","ft","s","C","°R","lb-mol","lm","rad","","")
dimtext(::typeof(normal(IPS))) = Values("lb","slinch","in","s","C","°R","slinch-mol","lm","rad","","")
dimtext(::typeof(normal(FPS))) = Values("pdl","lb","ft","s","C","°R","lb-mol","lm","rad","","")
dimtext(::typeof(normal(Gauss))) = Values("dyn","g","cm","s","","K","mol","lm","rad","","")

dimtext(::typeof(normal(IAU☉))) = Values("?","M☉","au","D","C","K","mol","lm","rad","","")
dimtext(::typeof(normal(IAUE))) = Values("?","ME","au","D","C","K","mol","lm","rad","","")
dimtext(::typeof(normal(IAUJ))) = Values("?","MJ","au","D","C","K","mol","lm","rad","","")
dimtext(::typeof(normal(MTS))) = Values("tf","t","m","s","C","K","mol","lm","rad","","")
dimtext(::typeof(normal(KKH))) = Values("kgf","kg","km","h","C","K","mol","lm","rad","","")
dimtext(::typeof(normal(MPH))) = Values("?","lb","mi","h","C","°R","lb-mol","lm","rad","","")
dimtext(::typeof(normal(Nautical))) = Values("kegf","keg","nm","h","eC","K","eg-mol","lm","rad","","")
dimtext(::typeof(normal(FFF))) = Values("?","fir","fur","ftn","Inf","°R","fir-mol","lm","rad","","")

dimtext(::typeof(normal(SI2019))) = dimtext(normal(Metric))
dimtext(::typeof(normal(SI1976))) = dimtext(normal(Metric))
dimtext(::typeof(normal(CODATA))) = dimtext(normal(Metric))
dimtext(::typeof(normal(Conventional))) = dimtext(normal(Metric))
dimtext(::typeof(normal(International))) = dimtext(normal(Metric))
dimtext(::typeof(normal(InternationalMean))) = dimtext(normal(Metric))
dimtext(::typeof(normal(SI2019Engineering))) = dimtext(normal(MetricEngineering))
dimtext(::typeof(normal(GravitationalSI2019))) = dimtext(normal(GravitationalMetric))
dimtext(::typeof(normal(MeridianEngineering))) = dimtext(normal(Meridian))
dimtext(::typeof(normal(GravitationalMeridian))) = dimtext(normal(Meridian))
#dimtext(::typeof(normal(British2019))) = dimtext(normal(British))
#dimtext(::typeof(normal(English2019))) = dimtext(normal(English))
#dimtext(::typeof(normal(IPS2019))) = dimtext(normal(IPS))
#dimtext(::typeof(normal(FPS2019))) = dimtext(normal(FPS))
dimtext(::typeof(normal(Survey))) = dimtext(normal(English))
#dimtext(::typeof(normal(Survey2019))) = dimtext(normal(English))
dimtext(::typeof(normal(EMU))) = dimtext(normal(Gauss))
dimtext(::typeof(normal(ESU))) = dimtext(normal(Gauss))
dimtext(::typeof(normal(LorentzHeaviside))) = dimtext(normal(Gauss))
#dimtext(::typeof(normal(Thomson))) = dimtext(normal(Gauss))

for U ∈ (:MetricEngineering, :SI2019Engineering,:GravitationalMetric,:GravitationalSI2019)
    @eval begin
showgroup(io::IO,::typeof($U(frequency)),::typeof(normal($U))) = print(io,"Hz")
showgroup(io::IO,::typeof($U(frequencydrift)),::typeof(normal($U))) = print(io,"Hz*s⁻¹")
showgroup(io::IO,::typeof($U(illuminance)),::typeof(normal($U))) = print(io,"lx")
showgroup(io::IO,::typeof(luminance),::typeof(normal($U))) = print(io,"nt")
showgroup(io::IO,::typeof($U(luminousexposure)),::typeof(normal($U))) = print(io,"lx*s")
    end
end
showgroup(io::IO,::typeof(MetricEngineering(specificforce)),::typeof(normal(MetricEngineering))) = print(io,"g₀")
showgroup(io::IO,::typeof(SI2019Engineering(specificforce)),::typeof(normal(SI2019Engineering))) = print(io,"g₀")
showgroup(io::IO,::typeof(MeridianEngineering(specificforce)),::typeof(normal(MeridianEngineering))) = print(io,"g₀")
showgroup(io::IO,::typeof(English(specificforce)),::typeof(normal(English))) = print(io,"g₀")
showgroup(io::IO,::typeof(Survey(specificforce)),::typeof(normal(Survey))) = print(io,"g₀")
for U ∈ (:Metric, :SI2019, :CODATA, :Conventional, :International, :InternationalMean)
    @eval begin
showgroup(io::IO,::typeof($U(frequency)),::typeof(normal($U))) = print(io,"Hz")
showgroup(io::IO,::typeof($U(frequencydrift)),::typeof(normal($U))) = print(io,"Hz*s⁻¹")
showgroup(io::IO,::typeof($U(force)),::typeof(normal($U))) = print(io,"N")
showgroup(io::IO,::typeof($U(inv(force))),::typeof(normal($U))) = print(io,"N⁻¹")
showgroup(io::IO,::typeof($U(pressure)),::typeof(normal($U))) = print(io,"Pa")
showgroup(io::IO,::typeof($U(energy)),::typeof(normal($U))) = print(io,"J")
showgroup(io::IO,::typeof($U(inv(energy))),::typeof(normal($U))) = print(io,"J⁻¹")
showgroup(io::IO,::typeof($U(power)),::typeof(normal($U))) = print(io,"W")
showgroup(io::IO,::typeof($U(inv(power))),::typeof(normal($U))) = print(io,"W⁻¹")
showgroup(io::IO,::typeof($U(electricpotential)),::typeof(normal($U))) = print(io,"V")
showgroup(io::IO,::typeof($U(inv(electricpotential))),::typeof(normal($U))) = print(io,"V⁻¹")
showgroup(io::IO,::typeof($U(capacitance)),::typeof(normal($U))) = print(io,"F")
showgroup(io::IO,::typeof($U(inv(capacitance))),::typeof(normal($U))) = print(io,"F⁻¹")
showgroup(io::IO,::typeof($U(resistance)),::typeof(normal($U))) = print(io,"Ω")
showgroup(io::IO,::typeof($U(conductance)),::typeof(normal($U))) = print(io,"S")
showgroup(io::IO,::typeof($U(magneticflux)),::typeof(normal($U))) = print(io,"Wb")
showgroup(io::IO,::typeof($U(magneticfluxdensity)),::typeof(normal($U))) = print(io,"T")
showgroup(io::IO,::typeof($U(permeance)),::typeof(normal($U))) = print(io,"H")
showgroup(io::IO,::typeof($U(luminousflux/power)),::typeof(normal($U))) = print(io,"lm*W⁻¹")
showgroup(io::IO,::typeof($U(luminousintensity)),::typeof(normal($U))) = print(io,"cd")
showgroup(io::IO,::typeof($U(illuminance)),::typeof(normal($U))) = print(io,"lx")
showgroup(io::IO,::typeof(luminance),::typeof(normal($U))) = print(io,"nt")
showgroup(io::IO,::typeof($U(luminousexposure)),::typeof(normal($U))) = print(io,"lx*s")
showgroup(io::IO,::typeof($U(catalysis)),::typeof(normal($U))) = print(io,"kat")

showgroup(io::IO,::typeof($U(impulse)),::typeof(normal($U))) = print(io,"N*s")
showgroup(io::IO,::typeof($U(angularmomentum)),::typeof(normal($U))) = print(io,"J*s")
showgroup(io::IO,::typeof($U(action*speed)),::typeof(normal($U))) = print(io,"J*m")
showgroup(io::IO,::typeof($U(yank)),::typeof(normal($U))) = print(io,"N*s⁻¹")
showgroup(io::IO,::typeof($U(fluence)),::typeof(normal($U))) = print(io,"N*m⁻¹")
showgroup(io::IO,::typeof($U(compressibility)),::typeof(normal($U))) = print(io,"Pa⁻¹")
showgroup(io::IO,::typeof($U(vectorpotential)),::typeof(normal($U))) = print(io,"Wb*m⁻¹")
showgroup(io::IO,::typeof($U(magneticmoment)),::typeof(normal($U))) = print(io,"Wb*m")
showgroup(io::IO,::typeof($U(mobility)),::typeof(normal($U))) = print(io,"m²s⁻¹V⁻¹")

showgroup(io::IO,::typeof($U(viscosity)),::typeof(normal($U))) = print(io,"Pa*s")
showgroup(io::IO,::typeof($U(intensity)),::typeof(normal($U))) = print(io,"W*m⁻²")
showgroup(io::IO,::typeof($U(powerdensity)),::typeof(normal($U))) = print(io,"W*m⁻³")
showgroup(io::IO,::typeof($U(intensity/Θ^4)),::typeof(normal($U))) = print(io,"W*m⁻²K⁻⁴")
showgroup(io::IO,::typeof($U(pressure/Θ^4)),::typeof(normal($U))) = print(io,"J*m⁻³K⁻⁴")
showgroup(io::IO,::typeof($U(𝟙/T/Θ)),::typeof(normal($U))) = print(io,"Hz*K⁻¹")
showgroup(io::IO,::typeof($U(entropy/Q)),::typeof(normal($U))) = print(io,"V*K⁻¹")
showgroup(io::IO,::typeof($U(entropy)),::typeof(normal($U))) = print(io,"J*K⁻¹")
showgroup(io::IO,::typeof($U(specificentropy)),::typeof(normal($U))) = print(io,"J*K⁻¹kg⁻¹")
showgroup(io::IO,::typeof($U(specificenergy)),::typeof(normal($U))) = print(io,"J*kg⁻¹")
showgroup(io::IO,::typeof($U(thermalconductivity)),::typeof(normal($U))) = print(io,"W*m⁻¹K⁻¹")
showgroup(io::IO,::typeof($U(thermalconductance)),::typeof(normal($U))) = print(io,"W*K⁻¹")
showgroup(io::IO,::typeof($U(thermalresistance)),::typeof(normal($U))) = print(io,"K*W⁻¹")
showgroup(io::IO,::typeof($U(thermalresistivity)),::typeof(normal($U))) = print(io,"K*m*W⁻¹")
showgroup(io::IO,::typeof($U(molarconductivity)),::typeof(normal($U))) = print(io,"S*m²mol⁻¹")
showgroup(io::IO,::typeof($U(electricpotential/M)),::typeof(normal($U))) = print(io,"V*kg⁻¹")
showgroup(io::IO,::typeof($U(inv(magneticflux))),::typeof(normal($U))) = print(io,"Hz*V⁻¹")
showgroup(io::IO,::typeof($U(action*speed/Q)),::typeof(normal($U))) = print(io,"V*m")
showgroup(io::IO,::typeof($U(electricfield)),::typeof(normal($U))) = print(io,"V*m⁻¹")
showgroup(io::IO,::typeof($U(permittivity)),::typeof(normal($U))) = print(io,"F*m⁻¹")
showgroup(io::IO,::typeof($U(inv(permittivity))),::typeof(normal($U))) = print(io,"m*F⁻¹")
showgroup(io::IO,::typeof($U(permeability)),::typeof(normal($U))) = print(io,"H*m⁻¹")
showgroup(io::IO,::typeof($U(inv(permeability))),::typeof(normal($U))) = print(io,"m*H⁻¹")
showgroup(io::IO,::typeof($U(resistivity)),::typeof(normal($U))) = print(io,"Ω*m")
showgroup(io::IO,::typeof($U(conductivity)),::typeof(normal($U))) = print(io,"S*m⁻¹")
showgroup(io::IO,::typeof($U(reluctance)),::typeof(normal($U))) = print(io,"H⁻¹")
showgroup(io::IO,::typeof($U(magneticdipolemoment)),::typeof(normal($U))) = print(io,"J*T⁻¹")
showgroup(io::IO,::typeof($U(molarenergy)),::typeof(normal($U))) = print(io,"J*mol⁻¹")
showgroup(io::IO,::typeof($U(molarentropy)),::typeof(normal($U))) = print(io,"J*K⁻¹mol⁻¹")
    end
end

showgroup(io::IO,::typeof(Meridian(frequency)),::typeof(normal(Meridian))) = print(io,"Hz")
showgroup(io::IO,::typeof(Meridian(frequencydrift)),::typeof(normal(Meridian))) = print(io,"Hz*s⁻¹")
showgroup(io::IO,::typeof(Meridian(force)),::typeof(normal(Meridian))) = print(io,"eN")
showgroup(io::IO,::typeof(Meridian(inv(force))),::typeof(normal(Meridian))) = print(io,"eN⁻¹")
showgroup(io::IO,::typeof(Meridian(pressure)),::typeof(normal(Meridian))) = print(io,"ePa")
showgroup(io::IO,::typeof(Meridian(energy)),::typeof(normal(Meridian))) = print(io,"eJ")
showgroup(io::IO,::typeof(Meridian(inv(energy))),::typeof(normal(Meridian))) = print(io,"eJ⁻¹")
showgroup(io::IO,::typeof(Meridian(power)),::typeof(normal(Meridian))) = print(io,"eW")
showgroup(io::IO,::typeof(Meridian(inv(power))),::typeof(normal(Meridian))) = print(io,"eW⁻¹")
showgroup(io::IO,::typeof(Meridian(electricpotential)),::typeof(normal(Meridian))) = print(io,"eV")
showgroup(io::IO,::typeof(Meridian(inv(electricpotential))),::typeof(normal(Meridian))) = print(io,"eV⁻¹")
showgroup(io::IO,::typeof(Meridian(capacitance)),::typeof(normal(Meridian))) = print(io,"eF")
showgroup(io::IO,::typeof(Meridian(inv(capacitance))),::typeof(normal(Meridian))) = print(io,"eF⁻¹")
showgroup(io::IO,::typeof(Meridian(resistance)),::typeof(normal(Meridian))) = print(io,"eΩ")
showgroup(io::IO,::typeof(Meridian(conductance)),::typeof(normal(Meridian))) = print(io,"eS")
showgroup(io::IO,::typeof(Meridian(magneticflux)),::typeof(normal(Meridian))) = print(io,"eWb")
showgroup(io::IO,::typeof(Meridian(magneticfluxdensity)),::typeof(normal(Meridian))) = print(io,"eT")
showgroup(io::IO,::typeof(Meridian(permeance)),::typeof(normal(Meridian))) = print(io,"eH")
showgroup(io::IO,::typeof(Meridian(luminousflux/power)),::typeof(normal(Meridian))) = print(io,"lm*eW⁻¹")
showgroup(io::IO,::typeof(Meridian(luminousintensity)),::typeof(normal(Meridian))) = print(io,"cd")
showgroup(io::IO,::typeof(Meridian(illuminance)),::typeof(normal(Meridian))) = print(io,"elx")
showgroup(io::IO,::typeof(luminance),::typeof(normal(Meridian))) = print(io,"ent")
showgroup(io::IO,::typeof(Meridian(luminousexposure)),::typeof(normal(Meridian))) = print(io,"lx*s")
showgroup(io::IO,::typeof(Meridian(catalysis)),::typeof(normal(Meridian))) = print(io,"ekat")

showgroup(io::IO,::typeof(Meridian(impulse)),::typeof(normal(Meridian))) = print(io,"eN*s")
showgroup(io::IO,::typeof(Meridian(angularmomentum)),::typeof(normal(Meridian))) = print(io,"eJ*s")
showgroup(io::IO,::typeof(Meridian(action*speed)),::typeof(normal(Meridian))) = print(io,"eJ*em")
showgroup(io::IO,::typeof(Meridian(yank)),::typeof(normal(Meridian))) = print(io,"eN*s⁻¹")
showgroup(io::IO,::typeof(Meridian(fluence)),::typeof(normal(Meridian))) = print(io,"eN*em⁻¹")
showgroup(io::IO,::typeof(Meridian(compressibility)),::typeof(normal(Meridian))) = print(io,"ePa⁻¹")
showgroup(io::IO,::typeof(Meridian(vectorpotential)),::typeof(normal(Meridian))) = print(io,"eWb*em⁻¹")
showgroup(io::IO,::typeof(Meridian(magneticmoment)),::typeof(normal(Meridian))) = print(io,"eWb*em")
showgroup(io::IO,::typeof(Meridian(mobility)),::typeof(normal(Meridian))) = print(io,"em²s⁻¹eV⁻¹")

showgroup(io::IO,::typeof(Meridian(viscosity)),::typeof(normal(Meridian))) = print(io,"ePa*s")
showgroup(io::IO,::typeof(Meridian(intensity)),::typeof(normal(Meridian))) = print(io,"eW*em⁻²")
showgroup(io::IO,::typeof(Meridian(powerdensity)),::typeof(normal(Meridian))) = print(io,"eW*m⁻³")
showgroup(io::IO,::typeof(Meridian(intensity/Θ^4)),::typeof(normal(Meridian))) = print(io,"eW*em⁻²K⁻⁴")
showgroup(io::IO,::typeof(Meridian(pressure/Θ^4)),::typeof(normal(Meridian))) = print(io,"eJ*em⁻³K⁻⁴")
showgroup(io::IO,::typeof(Meridian(𝟙/T/Θ)),::typeof(normal(Meridian))) = print(io,"Hz*K⁻¹")
showgroup(io::IO,::typeof(Meridian(entropy/Q)),::typeof(normal(Meridian))) = print(io,"eV*K⁻¹")
showgroup(io::IO,::typeof(Meridian(entropy)),::typeof(normal(Meridian))) = print(io,"eJ*K⁻¹")
showgroup(io::IO,::typeof(Meridian(specificentropy)),::typeof(normal(Meridian))) = print(io,"eJ*K⁻¹keg⁻¹")
showgroup(io::IO,::typeof(Meridian(specificenergy)),::typeof(normal(Meridian))) = print(io,"eJ*keg⁻¹")
showgroup(io::IO,::typeof(Meridian(thermalconductivity)),::typeof(normal(Meridian))) = print(io,"eW*em⁻¹K⁻¹")
showgroup(io::IO,::typeof(Meridian(thermalresistance)),::typeof(normal(Meridian))) = print(io,"K*eW⁻¹")
showgroup(io::IO,::typeof(Meridian(thermalresistivity)),::typeof(normal(Meridian))) = print(io,"K*em*eW⁻¹")
showgroup(io::IO,::typeof(Meridian(molarconductivity)),::typeof(normal(Meridian))) = print(io,"eS*em²eg-mol⁻¹")
showgroup(io::IO,::typeof(Meridian(electricpotential/M)),::typeof(normal(Meridian))) = print(io,"eV*kg⁻¹")
showgroup(io::IO,::typeof(Meridian(inv(magneticflux))),::typeof(normal(Meridian))) = print(io,"Hz*eV⁻¹")
showgroup(io::IO,::typeof(Meridian(action*speed/Q)),::typeof(normal(Meridian))) = print(io,"eV*em")
showgroup(io::IO,::typeof(Meridian(electricfield)),::typeof(normal(Meridian))) = print(io,"eV*em⁻¹")
showgroup(io::IO,::typeof(Meridian(permittivity)),::typeof(normal(Meridian))) = print(io,"eF*em⁻¹")
showgroup(io::IO,::typeof(Meridian(inv(permittivity))),::typeof(normal(Meridian))) = print(io,"em*eF⁻¹")
showgroup(io::IO,::typeof(Meridian(permeability)),::typeof(normal(Meridian))) = print(io,"eH*em⁻¹")
showgroup(io::IO,::typeof(Meridian(inv(permeability))),::typeof(normal(Meridian))) = print(io,"em*eH⁻¹")
showgroup(io::IO,::typeof(Meridian(resistivity)),::typeof(normal(Meridian))) = print(io,"eΩ*em")
showgroup(io::IO,::typeof(Meridian(conductivity)),::typeof(normal(Meridian))) = print(io,"eS*em⁻¹")
showgroup(io::IO,::typeof(Meridian(reluctance)),::typeof(normal(Meridian))) = print(io,"eH⁻¹")
showgroup(io::IO,::typeof(Meridian(magneticdipolemoment)),::typeof(normal(Meridian))) = print(io,"eJ*eT⁻¹")
showgroup(io::IO,::typeof(Meridian(molarenergy)),::typeof(normal(Meridian))) = print(io,"eJ*eg-mol⁻¹")
showgroup(io::IO,::typeof(Meridian(molarentropy)),::typeof(normal(Meridian))) = print(io,"eJ*K⁻¹eg-mol⁻¹")

for U ∈ (:Gauss, :EMU, :ESU, :LorentzHeaviside)
    @eval begin
showgroup(io::IO,::typeof($U(frequency)),::typeof(normal($U))) = print(io,"Hz")
showgroup(io::IO,::typeof($U(force)),::typeof(normal($U))) = print(io,"dyn")
showgroup(io::IO,::typeof($U(inv(force))),::typeof(normal($U))) = print(io,"dyn⁻¹")
showgroup(io::IO,::typeof($U(specificforce)),::typeof(normal($U))) = print(io,"gal")
showgroup(io::IO,::typeof($U(specificforce/L)),::typeof(normal($U))) = print(io,"gal*cm⁻¹")
showgroup(io::IO,::typeof($U(pressure)),::typeof(normal($U))) = print(io,"Ba")
showgroup(io::IO,::typeof($U(energy)),::typeof(normal($U))) = print(io,"erg")
showgroup(io::IO,::typeof($U(inv(energy))),::typeof(normal($U))) = print(io,"erg⁻¹")
showgroup(io::IO,::typeof($U(power)),::typeof(normal($U))) = print(io,"erg*s⁻¹")
showgroup(io::IO,::typeof($U(inv(power))),::typeof(normal($U))) = print(io,"s*erg⁻¹")

showgroup(io::IO,::typeof($U(luminousflux/power)),::typeof(normal($U))) = print(io,"lm*s*erg⁻¹")
showgroup(io::IO,::typeof($U(luminousintensity)),::typeof(normal($U))) = print(io,"cd")
showgroup(io::IO,::typeof($U(illuminance)),::typeof(normal($U))) = print(io,"ph")
showgroup(io::IO,::typeof(luminance),::typeof(normal($U))) = print(io,"sb")

showgroup(io::IO,::typeof($U(impulse)),::typeof(normal($U))) = print(io,"dyn*s")
showgroup(io::IO,::typeof($U(angularmomentum)),::typeof(normal($U))) = print(io,"erg*s")
showgroup(io::IO,::typeof($U(yank)),::typeof(normal($U))) = print(io,"dyn*s⁻¹")
showgroup(io::IO,::typeof($U(compressibility)),::typeof(normal($U))) = print(io,"Ba⁻¹")

showgroup(io::IO,::typeof($U(viscosity)),::typeof(normal($U))) = print(io,"P")
showgroup(io::IO,::typeof($U(diffusivity)),::typeof(normal($U))) = print(io,"St")
showgroup(io::IO,::typeof($U(pressure/Θ^4)),::typeof(normal($U))) = print(io,"Ba*K⁻⁴")
showgroup(io::IO,::typeof($U(𝟙/T/Θ)),::typeof(normal($U))) = print(io,"Hz*K⁻¹")
showgroup(io::IO,::typeof($U(entropy)),::typeof(normal($U))) = print(io,"erg*K⁻¹")
showgroup(io::IO,::typeof($U(specificentropy)),::typeof(normal($U))) = print(io,"erg*K⁻¹g⁻¹")
showgroup(io::IO,::typeof($U(specificenergy)),::typeof(normal($U))) = print(io,"erg*g⁻¹")
showgroup(io::IO,::typeof($U(thermalconductivity)),::typeof(normal($U))) = print(io,"erg*s⁻¹m⁻¹K⁻¹")
showgroup(io::IO,::typeof($U(thermalresistivity)),::typeof(normal($U))) = print(io,"K*m*s*erg⁻¹")
    end
end

showgroup(io::IO,::typeof(EMU(current)),::typeof(normal(EMU))) = print(io,"Bi")
showgroup(io::IO,::typeof(EMU(magneticflux)),::typeof(normal(EMU))) = print(io,"Mx")
showgroup(io::IO,::typeof(EMU(magneticfluxdensity)),::typeof(normal(EMU))) = print(io,"G")
showgroup(io::IO,::typeof(EMU(magneticfield)),::typeof(normal(EMU))) = print(io,"Oe")
showgroup(io::IO,::typeof(EMU(reluctance)),::typeof(normal(EMU))) = print(io,"Bi*Mx⁻¹")
showgroup(io::IO,::typeof(EMU(magneticdipolemoment)),::typeof(normal(EMU))) = print(io,"erg*G⁻¹")
showgroup(io::IO,::typeof(EMU(vectorpotential)),::typeof(normal(EMU))) = print(io,"Mx*cm⁻¹")
showgroup(io::IO,::typeof(EMU(magneticmoment)),::typeof(normal(EMU))) = print(io,"Mx*cm")

showgroup(io::IO,::typeof(Gauss(charge)),::typeof(normal(Gauss))) = print(io,"Fr")
showgroup(io::IO,::typeof(Gauss(magneticflux)),::typeof(normal(Gauss))) = print(io,"Mx")
showgroup(io::IO,::typeof(Gauss(magneticfluxdensity)),::typeof(normal(Gauss))) = print(io,"G")
showgroup(io::IO,::typeof(Gauss(magneticfield)),::typeof(normal(Gauss))) = print(io,"Oe")
showgroup(io::IO,::typeof(Gauss(reluctance)),::typeof(normal(Gauss))) = print(io,"Fr*s⁻¹Mx⁻¹")
showgroup(io::IO,::typeof(Gauss(magneticdipolemoment)),::typeof(normal(Gauss))) = print(io,"erg*G⁻¹")
showgroup(io::IO,::typeof(Gauss(vectorpotential)),::typeof(normal(Gauss))) = print(io,"Mx*cm⁻¹")
showgroup(io::IO,::typeof(Gauss(magneticmoment)),::typeof(normal(Gauss))) = print(io,"Mx*cm")

showgroup(io::IO,::typeof(GravitationalMetric(mass)),::typeof(normal(GravitationalMetric))) = print(io,"hyl")
showgroup(io::IO,::typeof(GravitationalSI2019(mass)),::typeof(normal(GravitationalSI2019))) = print(io,"hyl")
showgroup(io::IO,::typeof(GravitationalSI2019(mass)),::typeof(normal(GravitationalMeridian))) = print(io,"ehyl")
showgroup(io::IO,::typeof(British(mass)),::typeof(normal(British))) = print(io,"slug")
showgroup(io::IO,::typeof(IPS(mass)),::typeof(normal(IPS))) = print(io,"slinch")
showgroup(io::IO,::typeof(FPS(force)),::typeof(normal(FPS))) = print(io,"pdl")

showgroup(io::IO,::typeof(British(density)),::typeof(normal(British))) = print(io,"slug*ft⁻³")
showgroup(io::IO,::typeof(IPS(density)),::typeof(normal(IPS))) = print(io,"slinch*in⁻³")
showgroup(io::IO,::typeof(GravitationalMetric(density)),::typeof(normal(GravitationalMetric))) = print(io,"hyl*m⁻³")
showgroup(io::IO,::typeof(GravitationalSI2019(density)),::typeof(normal(GravitationalSI2019))) = print(io,"hyl*m⁻³")
showgroup(io::IO,::typeof(GravitationalMeridian(density)),::typeof(normal(GravitationalMeridian))) = print(io,"ehyl*m⁻³")

showgroup(io::IO,::typeof(English(luminousintensity)),::typeof(normal(English))) = print(io,"cd")
showgroup(io::IO,::typeof(English(illuminance)),::typeof(normal(English))) = print(io,"fc")
showgroup(io::IO,::typeof(British(luminousintensity)),::typeof(normal(British))) = print(io,"cd")
showgroup(io::IO,::typeof(British(illuminance)),::typeof(normal(British))) = print(io,"fc")
showgroup(io::IO,::typeof(FPS(luminousintensity)),::typeof(normal(FPS))) = print(io,"cd")
showgroup(io::IO,::typeof(FPS(illuminance)),::typeof(normal(FPS))) = print(io,"fc")

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
