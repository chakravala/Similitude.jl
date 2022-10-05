module Similitude

#   This file is part of Similitude.jl
#   It is licensed under the MIT license
#   Similitude Copyright (C) 2020 Michael Reed
#       _           _                         _
#      | |         | |                       | |
#   ___| |__   __ _| | ___ __ __ ___   ____ _| | __ _
#  / __| '_ \ / _` | |/ / '__/ _` \ \ / / _` | |/ _` |
# | (__| | | | (_| |   <| | | (_| |\ V / (_| | | (_| |
#  \___|_| |_|\__,_|_|\_\_|  \__,_| \_/ \__,_|_|\__,_|
#
#   https://github.com/chakravala
#   https://crucialflow.com
#     _     _         _  _  _                   _
#    | |   (_)       (_)| |(_) _               | |
#     \ \   _  ____   _ | | _ | |_   _   _   _ | |  ____
#      \ \ | ||    \ | || || ||  _) | | | | / || | / _  )
#  _____) )| || | | || || || || |__ | |_| |( (_| |( (/ /
# (______/ |_||_|_|_||_||_||_| \___) \____| \____| \____)

import Base: @pure
import UnitSystems
import UnitSystems: UnitSystem, universe, Coupling, logdb, expdb, dB
import UnitSystems: Systems, Dimensionless, Constants, Physics, Convert, Derived
const dir = dirname(pathof(UnitSystems))

using LinearAlgebra
import AbstractTensors: TupleVector, Variables, Values, value

import FieldConstants: param, Constant
import FieldAlgebra
import FieldAlgebra: showgroup, factorize, valueat, product, coefprod, makeint, measure
import FieldAlgebra: AbstractModule, AbelianGroup, Group, LogGroup, ExpGroup
import FieldAlgebra: value, isonezero, islog, base
#import FieldAlgebra: coef, coefprod, checkint, checkints, promoteint, promoteints, expos, chars,  makeint, findpower, latexpo, printexpo, printdims, printnum, showgroup, showfun, valueat

macro group(args...)
    FieldAlgebra.group(args...)
end

include("dimension.jl")
include("constant.jl")

dimtext(u) = isq
printdims(io::IO,x::Group{T,dims},u::UnitSystem) where T = printdims(io,x,dimtext(normal(u)))
printdims(io::IO,x::Group{T,vals},u::UnitSystem) where T = printdims(io,x,basis)
printdims(io::IO,x::Group{T,dims}) where T = printdims(io,x,isq)
printdims(io::IO,x::Group{T,vals}) where T = printdims(io,x,basis)

#logdb(x::Constant{D}) where D = Constant{logdb(D)}()
logdb(x::Quantity{D,U}) where {D,U} = Quantity{logdb(D),U}(logdb(x.v))

@pure UnitSystems.unit(x::Group{:USQ},y=1) = x
@pure UnitSystems.unit(x::AbelianGroup,y=1) = x
@pure UnitSystems.unit(x::ConvertUnit,y=1) = x

UnitSystems.unit(x::Group,y=1) = x
UnitSystems.unit(x::Quantity,y=1) = x

for unit ∈ (Systems...,Dimensionless...,Constants...,Physics...,Derived...,Convert...)
    unit ∉ (:length,:time,:angle) && @eval export $unit
end

@pure ratio_calc(D::ExpGroup{B},U,S) where B = B^ratio_calc(value(D),U,S)
@pure ratio_calc(D::LogGroup{B},U,S) where B = log(B,ratio_calc(value(D),U,S))
@pure ratio_calc(D::LogGroup{ℯ},U,S) = log(ratio_calc(value(D),U,S))
@pure ratio_calc(D::ExpGroup{ℯ},U,S) = exp(ratio_calc(value(D),U,S))
@pure ratio(::Constant{D},U,S) where D = ratio_calc(UnitSystem(D),normal(U),normal(S))
@pure function ratio_calc(D::Group,U,S)
    Constant(
    UnitSystems.boltzmann(U,S)^D.v[1]*
    UnitSystems.planckreduced(U,S)^D.v[2]*
    UnitSystems.lightspeed(U,S)^D.v[3]*
    UnitSystems.permeability(U,S)^D.v[4]*
    UnitSystems.electronmass(U,S)^D.v[5]*
    UnitSystems.molarmass(U,S)^D.v[6]*
    UnitSystems.luminousefficacy(U,S)^D.v[7]*
    UnitSystems.angle(U,S)^D.v[8]*
    UnitSystems.rationalization(U,S)^D.v[9]*
    UnitSystems.lorentz(U,S)^D.v[10]*
    UnitSystems.gravity(U,S)^D.v[11])
end
export Ratio
const Ratio = ratio

import UnitSystems: normal, unitname
normal(x::Quantity) = quantity(x)
(s::UnitSystem)(q::Quantity{D,U}) where {D,U} = (S=normal(s); Quantity{D,S}(q.v*ratio(D,U,S)))
(q::Quantity{D,U})(s::UnitSystem) where {D,U} = (S=normal(s); Quantity{D,S}(q.v*ratio(D,U,S)))

function Quantity(u::UnitSystem)
    U = constant(u); t = isone(normal(UnitSystems.boltzmann(U)))
    kB = Quantity{F*L/Θ,U}(UnitSystems.boltzmann(U))
    ħ = Quantity{F*L*T/A,U}(UnitSystems.planckreduced(U))
    c = Quantity{L/T,U}(UnitSystems.lightspeed(U))
    μ0 = Quantity{F*T*T*C*C/(Q*Q)/R,U}(UnitSystems.permeability(U))
    mₑ = Quantity{M,U}(UnitSystems.electronmass(U))
    Mᵤ = Quantity{M/N,U}(UnitSystems.molarmass(U))
    Kcd = Quantity{J*T/F/L,U}(UnitSystems.luminousefficacy(U))
    θ = Quantity{A,U}(UnitSystems.angle(U))
    λ = Quantity{R,U}(UnitSystems.rationalization(U))
    αL = Quantity{inv(C),U}(UnitSystems.lorentz(U))
    g₀ = Quantity{M*L/(F*T*T),U}(UnitSystems.gravity(U))
    τ = UnitSystems.tau(U)
    x = UnitSystems.two(U)
    y = UnitSystems.three(U)
    z = UnitSystems.five(U)
    u = UnitSystems.seven(U)
    v = UnitSystems.eleven(U)
    w = UnitSystems.nineteen(U)
    q = UnitSystems.fourtythree(U)
    UnitSystem(kB,ħ,c,μ0,mₑ,Mᵤ,Kcd,θ,λ,αL,g₀,Universe,τ,x,y,z,u,v,w,q)
end

const LD,JD = Constant(384399)*(𝟐*𝟓)^3,Constant(778479)*(𝟐*𝟓)^6
const μE☾ = Constant(UnitSystems.μE☾)

import UnitSystems: GaussSystem, EntropySystem, ElectricSystem, AstronomicalSystem
include("$dir/initdata.jl")

const Unified = Quantity(UnitSystem(F*L/Θ,F*L*T/A,L/T,F*T*T*C*C/(Q*Q)/R,M,M/N,J*T/F/L,A,R,inv(C),M*L/(F*T*T),Universe,τ,𝟐,𝟑,𝟓,𝟕,𝟏𝟏,𝟏𝟗,𝟒𝟑))
dimtext(::typeof(normal(Unified))) = Values("kB","ħ","𝘤","μ₀","mₑ","Mᵤ","Kcd","ϕ","λ","αL","g₀")

export Unified
unitname(::typeof(normal(Unified))) = "Unified"
(u::typeof(normal(Unified)))(d::Group) = UnitSystem(d)
(u::UnitSystem)(d::Group) = normal(Metric)(d)

for unit ∈ Convert
    if unit ∉ (:length,:time,:angle,:molarmass,:luminousefficacy)
        @eval const  $unit = dimensions(UnitSystems.$unit(UnitSystems.Natural,Natural))
    end
end
const gravityforce = acceleration/specificforce

@doc """
    (D::Constant)(U::UnitSystem,S::UnitSystem) = ConvertUnit{D,U,S}()

Constant for unit conversion for `D::Constant` from `U::UnitSystem` to `S::UnitSystem`.
```Julia
julia> energy(Metric,CGS)
$(energy(Metric,CGS))

julia> energy(Metric,English)
$(energy(Metric,English))
```
There still exists further opportunity to expand on the implementation of `ConvertUnit`.
""" ConvertUnit

@doc """
    (U::UnitSystem)(v::Number, D::Constant) ↦ Quantity(D,U,v) = Quantity{D,U}(v)

Numerical `Quantity` having value `v` with `D::Constant` specified in `U::UnitSystem`.
```Julia
julia> Metric(1,energy)
$(Metric(1,energy))

julia> English(1,energy)
$(English(1,energy))
```
An alternate syntax `Quantity(D::Constant, U::UnitSystem, v::Number)` is also available as standard syntax.
When `using UnitSystems` instead of `using Similitude`, this same syntax can be written so that code doesn't need to be changed while the output is generated.
""" Quantity

include("derived.jl")

(u::typeof(normal(UnitSystems.QCD)))(d::Group) = normal(Planck)(d)

for dim ∈ (:angle, :length, :time)
    @eval begin
        Base.sqrt(x::typeof($dim)) = sqrt(evaldim(x))
        Base.cbrt(x::typeof($dim)) = cbrt(evaldim(x))
        Base.inv(x::typeof($dim)) = inv(evaldim(x))
        Base.log(x::typeof($dim)) = log(evaldim(x))
        Base.log2(x::typeof($dim)) = log2(evaldim(x))
        Base.log10(x::typeof($dim)) = log10(evaldim(x))
        Base.log(b::Number,x::typeof($dim)) = log(b,evaldim(x))
        Base.exp(x::typeof($dim)) = exp(evaldim(x))
        Base.exp2(x::typeof($dim)) = exp2(evaldim(x))
        Base.exp10(x::typeof($dim)) = exp10(evaldim(x))
        Base.:^(b::T,x::typeof($dim)) where T<:Number = b^evaldim(x)
        Base.:^(a::typeof($dim),b::Number) = evaldim(a)^b
        Base.:^(a::typeof($dim),b::Integer) = evaldim(a)^b
        Base.:^(a::typeof($dim),b::Rational) = evaldim(a)^b
        Base.:*(a::typeof($dim),b::Constant) = evaldim(a)*b
        Base.:/(a::typeof($dim),b::Constant) = evaldim(a)/b
        Base.:*(a::Constant,b::typeof($dim)) = a*evaldim(b)
        Base.:/(a::Constant,b::typeof($dim)) = a/evaldim(b)
        Base.:*(a::Number,b::typeof($dim)) = a*evaldim(b)
        Base.:*(a::typeof($dim),b::Number) = evaldim(a)*b
        Base.:/(a::Number,b::typeof($dim)) = a*inv(evaldim(b))
        Base.:/(a::typeof($dim),b::Number) = evaldim(a)*inv(b)
    end
    for dim2 ∈ (:angle, :length, :time)
        @eval begin
            Base.:*(a::typeof($dim),b::typeof($dim2)) = evaldim(a)*evaldim(b)
            Base.:/(a::typeof($dim),b::typeof($dim2)) = evaldim(a)*evaldim(b)
        end
    end
end

if haskey(ENV,"UNITDOCS")
include("$dir/kinematicdocs.jl")
include("$dir/electromagneticdocs.jl")
include("$dir/thermodynamicdocs.jl")
include("$dir/physicsdocs.jl")
include("$dir/systems.jl")

@doc """
Physical dimension `Constant` represented by `Group` element `D`.
```Julia
F, M, L, T, Q, Θ, N, J, A, R, C
```
Operations on `Constant` are closed (`*`, `/`, `+`, `-`, `^`).
```Julia
julia> force(Unified)
$(force(Unified))

julia> mass(Unified)
$(mass(Unified))

julia> length(Unified)
$(length(Unified))

julia> time(Unified)
$(time(Unified))

julia> charge(Unified)
$(charge(Unified))

julia> temperature(Unified)
$(temperature(Unified))

julia> molaramount(Unified)
$(molaramount(Unified))

julia> luminousflux(Unified)
$(luminousflux(Unified))

julia> angle(Unified)
$(angle(Unified))

julia> rationalization(Unified)
$(rationalization(Unified))

julia> lorentz(Unified)
$(lorentz(Unified))
```
Derived dimension can be obtained from multiplicative base of 11 fundamental dimension symbols corresponding to `force`, `mass`, `length`, `time`, `charge`, `temperature`, `molaramount`, `luminousflux`, `angle`, `demagnetizingfactor`, and a `nonstandard` dimension.

Mechanics: `angle`, `$(listext(Kinematic))`, `$(listext(Mechanical))`;
Electromagnetics: `$(listext(Electromagnetic))`;
Thermodynamics: `$(listext(Thermodynamic))`,
`$(listext(Molar))`, `$(listext(Photometric))`.
""" USQ
export USQ

@doc """
    Unified = UnitSystem(...) # Unified System of Quantities (USQ)
    F, M, L, T, Q, Θ, N, J, A, R, C # fundamental base dimensions

Standard `Unified` system of `Quantities` (USQ) in terms of `UnitSystem` basis, transformed from the basis of `force`, `mass`, `length`, `time`, `charge`, `temperature`, `molaramount`, `luminousflux`, `angle`, `rationalization`, and a `nonstandard` dimension.

```Julia
julia> boltzmann(Unified) # entropy
$(boltzmann(Unified))

julia> planckreduced(Unified) # angularmomentum
$(planckreduced(Unified))

julia> lightspeed(Unified) # speed
$(lightspeed(Unified))

julia> vacuumpermeability(Unified) # permeability
$(vacuumpermeability(Unified))

julia> electronmass(Unified) # mass
$(electronmass(Unified))

julia> molarmass(Unified) # molarmass
$(molarmass(Unified))

julia> luminousefficacy(Unified) # luminousefficacy
$(luminousefficacy(Unified))

julia> radian(Unified) # angle
$(radian(Unified))

julia> rationalization(Unified) # demagnetizingfactor
$(rationalization(Unified))

julia> lorentz(Unified) # nonstandard
$(lorentz(Unified))

julia> gravity(Unified) # gravityforce
$(gravity(Unified))
```
""" Unified
end

end # module
