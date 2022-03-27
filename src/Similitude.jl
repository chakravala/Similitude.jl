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
import UnitSystems: UnitSystem, universe, Coupling, dB, Db
import UnitSystems: Systems, Dimensionless, Constants, Physics, Convert, Derived
const dir = dirname(pathof(UnitSystems))

using LinearAlgebra
import AbstractTensors: TupleVector, Values, value

include("group.jl")
include("dimension.jl")
include("constant.jl")

dimtext(u) = isq
printdims(io::IO,x::Group{T,dims},u::UnitSystem) where T = printdims(io,x,dimtext(normal(u)))
printdims(io::IO,x::Group{T,vals},u::UnitSystem) where T = printdims(io,x,basis)
printdims(io::IO,x::Group{T,dims}) where T = printdims(io,x,isq)
printdims(io::IO,x::Group{T,vals}) where T = printdims(io,x,basis)

dB(x::Dimension{D}) where D = Dimension{dB(D)}()
dB(x::Quantity{D,U}) where {D,U} = Quantity{dB(D),U}(dB(x.v))

@pure UnitSystems.unit(x::Dimension,y=1) = x
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
@pure ratio(::Dimension{D},U,S) where D = ratio_calc(UnitSystem(D),normal(U),normal(S))
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
    U = Constant(u); t = isone(normal(UnitSystems.boltzmann(U)))
    kB = Quantity{F*L/Θ,U}(UnitSystems.boltzmann(U))
    ħ = Quantity{F*L*T/A,U}(UnitSystems.planckreduced(U))
    c = Quantity{L/T,U}(UnitSystems.lightspeed(U))
    μ0 = Quantity{F*T*T*C*C/(Q*Q*A*A)/Λ,U}(UnitSystems.permeability(U))
    mₑ = Quantity{M,U}(UnitSystems.electronmass(U))
    Mᵤ = Quantity{M/N,U}(UnitSystems.molarmass(U))
    Kcd = Quantity{J*T/F/L,U}(UnitSystems.luminousefficacy(U))
    a = Quantity{𝟙,U}(UnitSystems.angle(U))
    λ = Quantity{Λ*A*A,U}(UnitSystems.rationalization(U))
    αL = Quantity{inv(C),U}(UnitSystems.lorentz(U))
    g₀ = Quantity{M*L/(F*T*T),U}(UnitSystems.gravity(U))
    τ = Quantity{A,U}(UnitSystems.twopi(U))
    x = UnitSystems.two(U)
    y = UnitSystems.three(U)
    z = UnitSystems.five(U)
    u = UnitSystems.seven(U)
    v = UnitSystems.eleven(U)
    w = UnitSystems.nineteen(U)
    q = UnitSystems.fourtythree(U)
    UnitSystem(kB,ħ,c,μ0,mₑ,Mᵤ,Kcd,a,λ,αL,g₀,Universe,τ,x,y,z,u,v,w,q)
end

const LD = Constant(UnitSystems.LD)
const μE☾ = Constant(UnitSystems.μE☾)

import UnitSystems: GaussSystem, EntropySystem, ElectricSystem, AstronomicalSystem
include("$dir/initdata.jl")

(u::UnitSystem)(d::Group) = normal(Metric)(d)

for unit ∈ Convert
    if unit ∉ (:length,:time,:angle,:molarmass,:luminousefficacy,:angularfrequency,:angularwavenumber,:angularmomentum,:solidangle)
        @eval const  $unit = dimensions(UnitSystems.$unit(UnitSystems.Natural,Natural))
    end
end
const solidangle = A*A
const angularfrequency = A/T
const angularwavenumber = A/L
const angularmomentum = F*L*T*A

include("derived.jl")

(u::typeof(normal(UnitSystems.QCD)))(d::Group) = normal(Planck)(d)

include("$dir/kinematicdocs.jl")
include("$dir/electromagneticdocs.jl")
include("$dir/thermodynamicdocs.jl")
include("$dir/physicsdocs.jl")
include("$dir/systems.jl")

end # module
