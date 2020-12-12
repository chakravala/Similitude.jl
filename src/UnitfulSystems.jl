module UnitfulSystems

#   This file is part of UnitfulSystems.jl. It is licensed under the MIT license
#   UnitfulSystems Copyright (C) 2020 Michael Reed

import Base: @pure
using UnitSystems, Unitful
import UnitSystems: Systems, Constants, Physics, Convert
export UnitSystems, Unitful, unitful
const ftlb = u"slug*ft^2/s^2"

UnitSystems.unit(x::Quantity,y=1) = x

for unit ∈ (Systems...,Constants...,Physics...,Convert...)
    unit ∉ (:length,:time) && @eval export $unit
end

for unit ∈ (Constants...,Physics...)
    if unit ∈ (:molarmass,:permeability,:permittivity,:charge,:magneticflux,:impedance,:conductance,:luminousefficacy)
        @eval @pure $unit(U::UnitSystem) = UnitSystems.$unit(U)
        @eval @doc $(string(@eval @doc UnitSystems.$unit)) $unit
    else
        @eval import UnitSystems.$unit
    end
end

for unit ∈ Convert
    @eval begin
        @pure @inline $unit(v::Number,U::UnitSystem) = $unit(v,U,Metric)
        @pure @inline $unit(v::Number,U::UnitSystem,S::UnitSystem) = ustrip((u=ustrip($unit(U,S));isone(u) ? v : v/u))*unit($unit(Natural,U))
        @pure @inline $unit(v::Number,U::UnitSystem{kB,ħ,𝘤,μ₀,mₑ},S::UnitSystem{kB,ħ,𝘤,μ₀,mₑ}) where {kB,ħ,𝘤,μ₀,mₑ} = ustrip(v)*unit($unit(UnitSystems.Natural,U))
        @pure @inline $unit(U::UnitSystem,S::UnitSystem) = UnitSystems.$unit(U,S)
    end
    if unit ∉ (Constants...,:permittivity,:charge,:magneticflux,:impedance,:conductance)
        @eval @pure @inline $unit(U::UnitSystem) = $unit(U,Metric)
        @eval @doc $(string(@eval @doc UnitSystems.$unit)) $unit
    end
end
for unit ∈ (:(Base.length),:(Base.time))
    @eval begin
        @pure @inline $unit(v::Quantity,U::UnitSystem) = $unit(v,U,Metric)
        @pure @inline $unit(v::Quantity,U::UnitSystem,S::UnitSystem) = ustrip((u=ustrip($unit(U,S));isone(u) ? v : v/u))*unit($unit(Natural,U))
        @pure @inline $unit(v::Quantity,U::UnitSystem{kB,ħ,𝘤,μ₀,mₑ},S::UnitSystem{kB,ħ,𝘤,μ₀,mₑ}) where {kB,ħ,𝘤,μ₀,mₑ} = ustrip(v)*unit($unit(Natural,U))
    end
end

"""
    unitful(::UnitSystem,JK=u"J/K",Js=u"J*s",ms=u"m/s",Hm=u"H/m",kg=u"kg")

Convert a `UnitSystem` to use `Unitful` values instead of plain numerical values.
"""
unitful(U::UnitSystem,JK=u"J/K",Js=u"J*s",ms=u"m/s",Hm=u"H/m",kg=u"kg") = U(JK,Js,ms,Hm,kg)

const EMU2019 = unitful(UnitSystems.EMU2019,u"erg/K",u"erg*s",u"cm/s",u"nH/cm",u"g")
const MTS = unitful(UnitSystems.MTS,u"kJ/K",u"kJ*s",u"m/s",u"kH/m",u"Mg")
const Metric = unitful(UnitSystems.Metric)
const SI2019 = unitful(UnitSystems.SI2019)
const CODATA = unitful(UnitSystems.CODATA)
const Conventional = unitful(UnitSystems.Conventional)
const English = unitful(UnitSystems.English,ftlb/u"Ra",ftlb*u"s",u"ft/s",1,u"slug")
const EnglishUS = unitful(UnitSystems.EnglishUS,ftlb/u"Ra",ftlb*u"s",u"ft/s",1,u"slug")

const SI = SI2019

@pure molarmass(U::UnitSystem{boltzmann(MTS)}) = molarmass(UnitSystems.MTS)*u"Mg/mol"
@pure molarmass(U::UnitSystem{boltzmann(EMU2019)}) = molarmass(UnitSystems.EMU2019)*u"g/mol"
@pure molarmass(U::UnitSystem{boltzmann(Metric)}) = molarmass(UnitSystems.Metric)*u"kg/mol"
@pure molarmass(U::UnitSystem{boltzmann(SI2019)}) = molarmass(UnitSystems.SI2019)*u"kg/mol"
@pure molarmass(U::UnitSystem{boltzmann(CODATA)}) = molarmass(UnitSystems.CODATA)*u"kg/mol"
@pure molarmass(U::UnitSystem{boltzmann(Conventional)}) = molarmass(UnitSystems.Conventional)*u"kg/mol"
@pure molarmass(U::UnitSystem{boltzmann(English)}) = molarmass(UnitSystems.English)*u"slug/mol"
@pure molarmass(U::UnitSystem{boltzmann(EnglishUS)}) = molarmass(UnitSystems.EnglishUS)*u"slug/mol"

for us ∈ (:EMU2019,:MTS,:Metric,:SI2019,:CODATA,:Conventional,:English,:EnglishUS)
    @eval @pure hyperfine(U::UnitSystem{boltzmann($us)}) = hyperfine(UnitSystems.$us)*u"Hz"
end
for us ∈ (:Metric,:SI2019,:CODATA,:Conventional)
    @eval @pure luminousefficacy(U::UnitSystem{boltzmann($us)}) = luminousefficacy(UnitSystems.$us)*u"cd/W"
end
@pure luminousefficacy(U::UnitSystem{boltzmann(EMU2019)}) = luminousefficacy(UnitSystems.EMU2019)*u"cd*s/erg"
@pure luminousefficacy(U::UnitSystem{boltzmann(MTS)}) = luminousefficacy(UnitSystems.MTS)*u"cd/GW"
@pure luminousefficacy(U::UnitSystem{boltzmann(English)}) = luminousefficacy(UnitSystems.English)*u"cd*s^3/slug/ft^2"
@pure luminousefficacy(U::UnitSystem{boltzmann(EnglishUS)}) = luminousefficacy(UnitSystems.EnglishUS)*u"cd*s^3/slug/ft^2"

for CAL ∈ (:cal,:calₜₕ,:cal₄,:cal₁₀,:cal₂₀,:calₘ,:calᵢₜ)
    KCAL = Symbol(:k,CAL)
    @eval const $CAL = UnitSystems.$CAL*u"cal"
    @eval const $KCAL = UnitSystems.$KCAL*u"kcal"
end

const atm = UnitSystems.atm*u"kPa"
const g₀ = UnitSystems.g₀*u"m/s^2"
const lbm = UnitSystems.lbm*u"ft/s^2"
const slug = UnitSystems.slug*u"kg/slug"
const ft = UnitSystems.ft*u"m/ft"
const ftUS = UnitSystems.ftUS*u"m/ft"
const rankine = UnitSystems.rankine*u"K/Ra"
const kelvin = UnitSystems.kelvin*u"Ra/K"

const ΔνCs = UnitSystems.ΔνCs*u"Hz"
const Kcd = UnitSystems.Kcd*u"cd/W"
const mP = UnitSystems.mP*u"kg"
const NA = UnitSystems.NA*u"mol^-1"
const kB = UnitSystems.kB*u"J/K"
const 𝘩 = UnitSystems.𝘩*u"J*s"
const 𝘤 = UnitSystems.𝘤*u"m/s"
const 𝘦 = UnitSystems.𝘦*u"C"
const R∞ = UnitSystems.R∞*u"m^-1"
const μ₀ = UnitSystems.μ₀*u"H/m"
const ħ = UnitSystems.ħ*u"J*s" # u"J*s/rad" ?
const δμ₀ = UnitSystems.δμ₀*u"H/m"
const Rᵤ = UnitSystems.Rᵤ*u"J/K/mol"
const mₑ = UnitSystems.mₑ*u"kg"
const RK1990 = UnitSystems.RK1990*u"Ω"
const RK2014 = UnitSystems.RK2014*u"Ω"
const KJ1990 = UnitSystems.KJ1990*u"Hz/V"
const KJ2014 = UnitSystems.KJ2014*u"Hz/V"
const ħ1990 = UnitSystems.ħ1990*u"J*s"
const ħ2014 = UnitSystems.ħ2014*u"J*s"
const mₑ1990 = UnitSystems.mₑ1990*u"kg"
const mₑ2014 = UnitSystems.mₑ2014*u"kg"

const GG = UnitSystems.GG*u"m^3/kg/s^2"
const κ = UnitSystems.κ*u"s^2/m/kg"
const σ = UnitSystems.σ*u"W/m^2/K^4" #
const μB = UnitSystems.μB*u"J/T" #
const ε₀ = UnitSystems.ε₀*u"F/m" #
const kₑ = UnitSystems.kₑ*u"N*m^2/C^2" #
const mₚ = UnitSystems.mₚ*u"kg"
const mᵤ = UnitSystems.mᵤ*u"kg"
const Mᵤ = UnitSystems.Mᵤ*u"kg/mol" #
const 𝔉 = UnitSystems.𝔉*u"C/mol" #
const Φ₀ = UnitSystems.Φ₀*u"Wb" #
const Z₀ = UnitSystems.Z₀*u"Ω" #
const G₀ = UnitSystems.G₀*u"S" #
const Eₕ = UnitSystems.Eₕ*u"J"
const a₀ = UnitSystems.a₀*u"m"
const rₑ = UnitSystems.rₑ*u"m"
const RK = UnitSystems.RK*u"Ω" #
const KJ = UnitSystems.KJ*u"Hz/V" #
const RH = UnitSystems.RH*u"1/m"
const Ry = UnitSystems.Ry*u"J"

const ℓP = UnitSystems.ℓP*u"m"
const tP = UnitSystems.tP*u"s"
const TP = UnitSystems.TP*u"K"

const lS = UnitSystems.lS*u"m"
const tS = UnitSystems.tS*u"s"
const mS = UnitSystems.mS*u"kg"
const qS = UnitSystems.qS*u"C"

const lA = UnitSystems.lA*u"m"
const tA = UnitSystems.tA*u"s"
const mA = UnitSystems.mA*u"kg"
const qA = UnitSystems.qA*u"C"

const lQCD = UnitSystems.lQCD*u"m"
const tQCD = UnitSystems.tQCD*u"s"
const mQCD = UnitSystems.mQCD*u"kg"

const Mu,Ru,SB,hh,cc,m0,e0,ke,me,mp,mu,ee,FF,Z0,G0,Eh,a0,re,g0,lP,ϵ₀,mB = Mᵤ,Rᵤ,σ,𝘩,𝘤,μ₀,ε₀,kₑ,mₑ,mₚ,mᵤ,𝘦,𝔉,Z₀,G₀,Eₕ,a₀,rₑ,g₀,ℓP,ε₀,μB
export κ, GG, NA, kB, Rᵤ, σ, 𝘩, ħ, 𝘤, μ₀, ε₀, kₑ, mₑ, mₚ, mᵤ, 𝘦, 𝔉, Φ₀, Z₀, G₀, Eₕ, R∞, a₀, rₑ, KJ, RK, Ru, SB, hh, cc, m0, e0, ke, me, mp, mu, ee, FF, Z0, G0, Eh, a0, re, μB
export αG, αinv, μₚₑ, μₑᵤ, μₚᵤ, mpe, meu, mpu, mP, δμ₀, Mᵤ, Mu, RH, Ry, ΔνCs, Kcd, ainv
export cal, kcal, calₜₕ, kcalₜₕ, calᵢₜ, kcalᵢₜ, ℓP, g₀, g0, atm, lbm, aG, BTUJ, BTUftlb
export lP, tP, TP, lS, tS, mS, qS, lA, tA, mA, qA, lQCD, tQCD, mQCD, ϵ₀, aL, αL

export slug, ft, KJ1990, KJ2014, RK1990, RK2014, mₑ1990, mₑ2014, temp, units
export slugs, kilograms, lbm, meters, feet, rankine, kelvin, moles, molecules
export UnitSystem, US, SI, CGS, CGS2019, CGSm, CGSe, HLU, FFF

@doc """
    EMU2019::UnitSystem{1e7*kB,1e7*ħ,100𝘤,1e7*μ₀,1000mₑ}

Centimetre-gram-second `UnitSystem` variant of tuned `SI2019` based on EMU (rationalized).

```Julia
julia> boltzmann(EMU2019)
$(boltzmann(EMU2019))

julia> planckreduced(EMU2019)
$(planckreduced(EMU2019))

julia> lightspeed(EMU2019)
$(lightspeed(EMU2019))

julia> permeability(EMU2019)
$(permeability(EMU2019))

julia> electronmass(EMU2019)
$(electronmass(EMU2019))
```
""" EMU2019

@doc """
    MTS::UnitSystem{1e6*Rᵤ*mₑ/μₑᵤ,1000ħ,𝘤,4π/1e4,mₑ/1000}

Metre-tonne-second `UnitSystem` variant of `Metric` system.

```Julia
julia> boltzmann(MTS)
$(boltzmann(MTS))

julia> planckreduced(MTS)
$(planckreduced(MTS))

julia> lightspeed(MTS)
$(lightspeed(MTS))

julia> permeability(MTS)
$(permeability(MTS))

julia> electronmass(MTS)
$(electronmass(MTS))
```
""" MTS

@doc """
    Metric::UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001,ħ,𝘤,4π*1e-7,mₑ}

Systeme International d'Unites (the SI units) adopted as the preferred `UnitSystem`.

```Julia
julia> boltzmann(Metric)
$(boltzmann(Metric))

julia> planckreduced(Metric)
$(planckreduced(Metric))

julia> lightspeed(Metric)
$(lightspeed(Metric))

julia> permeability(Metric)
$(permeability(Metric))

julia> electronmass(Metric)
$(electronmass(Metric))
```
""" Metric

@doc """
    SI2019::UnitSystem{kB,ħ,𝘤,μ₀,mₑ}

Systeme International d'Unites (the SI units) with `μ₀` for a tuned `charge` exactly.

```Julia
julia> boltzmann(SI2019)
$(boltzmann(SI2019))

julia> planckreduced(SI2019)
$(planckreduced(SI2019))

julia> lightspeed(SI2019)
$(lightspeed(SI2019))

julia> permeability(SI2019)
$(permeability(SI2019))

julia> electronmass(SI2019)
$(electronmass(SI2019))
```
""" SI2019, SI

@doc """
    CODATA::UnitSystem{Rᵤ*mₑ2014/μₑᵤ/0.001,2/RK2014/KJ2014^2/π,𝘤,2RK2014/𝘤/αinv,mₑ2014}

Metric `UnitSystem` based on Committee on Data of the International Science Council.

```Julia
julia> boltzmann(CODATA)
$(boltzmann(CODATA))

julia> planckreduced(CODATA)
$(planckreduced(CODATA))

julia> lightspeed(CODATA)
$(lightspeed(CODATA))

julia> permeability(CODATA)
$(permeability(CODATA))

julia> electronmass(CODATA)
$(electronmass(CODATA))
```
""" CODATA

@doc """
    Conventional::UnitSystem{Rᵤ*mₑ1990/μₑᵤ/0.001,2/RK1990/KJ1990^2/π,𝘤,2RK1990/𝘤/αinv,mₑ1990}

Conventional electronic `UnitSystem` with 1990 tuned `josephson` and `klitzing` constants.

```Julia
julia> boltzmann(Conventional)
$(boltzmann(Conventional))

julia> planckreduced(Conventional)
$(planckreduced(Conventional))

julia> lightspeed(Conventional)
$(lightspeed(Conventional))

julia> permeability(Conventional)
$(permeability(Conventional))

julia> electronmass(Conventional)
$(electronmass(Conventional))
```
""" Conventional

@doc """
    English::UnitSystem{kB*rankine/slug/ft^2,ħ/slug/ft^2,𝘤/ft,4π,mₑ/slug}

Engineering `UnitSystem` historically used by Britain and United States.

```Julia
julia> boltzmann(English)
$(boltzmann(English))

julia> planckreduced(English)
$(planckreduced(English))

julia> lightspeed(English)
$(lightspeed(English))

julia> electronmass(English)
$(electronmass(English))
```
""" English

@doc """
    EnglishUS::UnitSystem{1000Rᵤ*mₑ/μₑᵤ*rankine/slug/ftUS^2,ħ/slug/ftUS^2,𝘤/ftUS,4π,mₑ/slug}

Engineering `UnitSystem` based on the geophysical US survey foot (1200/3937).

```Julia
julia> boltzmann(EnglishUS)
$(boltzmann(EnglishUS))

julia> planckreduced(EnglishUS)
$(planckreduced(EnglishUS))

julia> lightspeed(EnglishUS)
$(lightspeed(EnglishUS))

julia> electronmass(EnglishUS)
$(electronmass(EnglishUS))
```
""" EnglishUS

end # module
