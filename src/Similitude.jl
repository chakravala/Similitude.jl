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
import UnitSystems: UnitSystem, universe, Coupling, logdb, expdb, dB, cache
import UnitSystems: Systems, Dimensionless, Constants, Physics, Convert, Derived
const dir = dirname(pathof(UnitSystems))

using LinearAlgebra

import FieldConstants: param, Constant
import FieldAlgebra
import FieldAlgebra: TupleVector, Variables, Values, value
import FieldAlgebra: showgroup, factorize, valueat, product, coefprod, makeint, measure
import FieldAlgebra: AbstractModule, AbelianGroup, Group, LogGroup, ExpGroup
import FieldAlgebra: value, isonezero, islog, base
#import FieldAlgebra: coef, coefprod, checkint, checkints, promoteint, promoteints, expos, chars,  makeint, findpower, latexpo, printexpo, printdims, printnum, showgroup, showfun, valueat

const CONSTDIM,CONSTVAL = false,false

macro group(args...)
    CONSTDIM ? FieldAlgebra.group(args...) : FieldAlgebra.group2(args...)
end
macro group2(args...)
    CONSTVAL ? FieldAlgebra.group(args...) : FieldAlgebra.group2(args...)
end

include("dimension.jl")
include("constant.jl")

dimtext(u) = isq
printdims(io::IO,x::Group{T,dims},u::UnitSystem) where T = printdims(io,x.v,dimtext(normal(u)))
printdims(io::IO,x::Group{T,vals},u::UnitSystem) where T = printdims(io,x.v,basis)
printdims(io::IO,x::Group{T,dims}) where T = printdims(io,x.v,isq)
printdims(io::IO,x::Group{T,vals}) where T = printdims(io,x.v,basis)

#logdb(x::Constant{D}) where D = Constant{logdb(D)}()
logdb(x::Quantity{U}) where U = Quantity{U}(logdb(x.v),logdb(dimensions(x)))

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
@pure ratio(D,U,S) = ratio_calc(UnitSystem(D),normal(U),normal(S))
@pure function ratio_calc(D::Group,U,S)
    (CONSTVAL ? Constant : identity)(
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
export Ratio, dimensions
const Ratio = ratio

import UnitSystems: normal, unitname
normal(x::Quantity) = quantity(x)
function (s::UnitSystem)(q::Quantity{U}) where U
    S = normal(s); D = dimensions(q)
    Quantity{S}(q.v*ratio(D,U,S),D)
end
function (q::Quantity{U})(s::UnitSystem) where U
    S = normal(s); D = dimensions(q)
    Quantity{S}(q.v*ratio(D,U,S),D)
end

function Quantity(u::UnitSystem)
    U = constant(u)
    kB = Quantity{U}(UnitSystems.boltzmann(U),F*L/Θ)
    ħ = Quantity{U}(UnitSystems.planckreduced(U),F*L*T/A)
    c = Quantity{U}(UnitSystems.lightspeed(U),L/T)
    μ0 = Quantity{U}(UnitSystems.permeability(U),F*T*T*C*C/(Q*Q)/R)
    mₑ = Quantity{U}(UnitSystems.electronmass(U),M)
    Mᵤ = Quantity{U}(UnitSystems.molarmass(U),M/N)
    Kcd = Quantity{U}(UnitSystems.luminousefficacy(U),J*T/F/L)
    θ = Quantity{U}(UnitSystems.angle(U),A)
    λ = Quantity{U}(UnitSystems.rationalization(U),R)
    αL = Quantity{U}(UnitSystems.lorentz(U),inv(C))
    g₀ = Quantity{U}(UnitSystems.gravity(U),M*L/(F*T*T))
    τ = UnitSystems.tau(U)
    x = UnitSystems.two(U)
    y = UnitSystems.three(U)
    z = UnitSystems.five(U)
    u = UnitSystems.seven(U)
    v = UnitSystems.eleven(U)
    w = UnitSystems.nineteen(U)
    q = UnitSystems.fourtythree(U)
    unitsystem(kB,ħ,c,μ0,mₑ,Mᵤ,Kcd,θ,λ,αL,g₀,Universe,τ,x,y,z,u,v,w,q)
end

const LD,JD = Constant(384399)*(𝟐*𝟓)^3,Constant(778479)*(𝟐*𝟓)^6
const μE☾ = Constant(UnitSystems.μE☾)

import UnitSystems: GaussSystem, EntropySystem, ElectricSystem, AstronomicalSystem

function includereplace(mod,file,str="Constant("=>"identity(")
    if VERSION >= v"1.6"
        for expr in Meta.parseall(replace(read(file,String),str)).args
            Base.eval(mod, expr)
        end
    else
        code = replace(read(file,String),str)
        pos = 1
        while pos <= lastindex(code)
            expr, next_pos = Meta.parse(code,pos; raise=false)
            expr === nothing && break # no more expressions
            Base.eval(mod,expr)
            pos = next_pos
        end
    end
end

println("Similitude: initializing UnitSystems data")
if CONSTVAL
    include("$dir/initdata.jl")
else
    includereplace(Similitude,"$dir/initdata.jl")
end

const Unified = Quantity(UnitSystem(F*L/Θ,F*L*T/A,L/T,F*T*T*C*C/(Q*Q)/R,M,M/N,J*T/F/L,A,R,inv(C),M*L/(F*T*T),Universe,τ,𝟐,𝟑,𝟓,𝟕,𝟏𝟏,𝟏𝟗,𝟒𝟑))
dimtext(::typeof(normal(Unified))) = Values("kB","ħ","𝘤","μ₀","mₑ","Mᵤ","Kcd","ϕ","λ","αL","g₀")
dimlatex(::typeof(normal(Unified))) = Values("\\text{k}_\\text{B}","\\hbar","\\text{c}","\\mu_0","\\text{m}_\\text{e}","\\text{M}_\\text{u}","\\text{K}_\\text{cd}","\\phi","\\lambda","\\alpha_\\text{L}","\\text{g}_0")

export Unified
unitname(::typeof(normal(Unified))) = "Unified"
(u::typeof(normal(Unified)))(d::Group) = d#UnitSystem(d)

for unit ∈ Convert
    if unit ∉ (:dimensionless,:length,:time,:angle,:molarmass,:luminousefficacy)
        @eval const  $unit = dimensions(UnitSystems.$unit(UnitSystems.Natural,Natural))
    end
end
const gravityforce = acceleration/specificforce

@doc """
    (D::Group{:USQ})(U::UnitSystem,S::UnitSystem) = ConvertUnit{U,S}(D)

Constant for unit conversion for `D::Group{:USQ}` from `U::UnitSystem` to `S::UnitSystem`.
```Julia
julia> energy(Metric,CGS)
$(energy(Metric,CGS))

julia> energy(Metric,English)
$(energy(Metric,English))
```
There still exists further opportunity to expand on the implementation of `ConvertUnit`.
""" ConvertUnit

@doc """
    (U::UnitSystem)(v::Number, d) ↦ Quantity(U,v,d) = Quantity{U}(v,d)

Numerical `Quantity` having value `v` with dimension `d` specified in `U::UnitSystem`.
```Julia
julia> Metric(1,energy)
$(Metric(1,energy))

julia> English(1,energy)
$(English(1,energy))
```
An alternate syntax `Quantity(U::UnitSystem, v::Number, d)` is also available as standard syntax.
When `using UnitSystems` instead of `using Similitude`, this same syntax can be written so that code doesn't need to be changed while the output is generated.
""" Quantity

dimlatex(U) = Values("\\text{F}","\\text{M}","\\text{L}","\\text{T}","\\text{Q}","\\Theta","\\text{N}","\\text{J}","\\text{A}","\\text{R}","\\text{C}")

println("Similitude: deriving Quantity measurements")
include("derived.jl")
include("unitdim.jl")

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

param1(x::Constant) = param(x)
param1(x) = x
compare(a::Symbol,b::Symbol,U) = param1(evaldim(a,U)) == param1(evaldim(b,U))
Base.:/(U::UnitSystem,::typeof(~)) = quotient(U)

export quotient
function quotient(U)
    len = length(Convert)
    out = [:x=>Symbol[] for i ∈ 1:len]
    for i ∈ 1:len
        x = Convert[i]
        out[i] = x=>Convert[findall(compare.(Ref(x),Convert,Ref(U)))]
    end
    i,num = 1,length(out)
    while i ≤ num
        for x ∈ out[i].second[2:end]
            j = findfirst(z->z==x,getproperty.(out,:first))
            if !isnothing(j) && j>i
                deleteat!(out,j)
                num -= 1
            end
        end
        i += 1
    end
    [evaldim(x.first,U)=>x.second for x ∈ out]
end

function printquotient(U)
    for x ∈ quotient(U)
        print("    $(x.first) => ")
        for j ∈ 1:length(x.second)
            y = x.second[j]
            print("$y ($(evaldim(y)))")
            j<length(x.second) && print(", ")
        end
        println()
    end
end

function latexquotient(U)
    latexquotient.(quotient(U),Ref(normal(U)))
end
function latexquotient(x,U)
    io = IOBuffer()
    Similitude.latexgroup(io,x.first,U)
    str1 = String(take!(io))
    for j ∈ 1:length(x.second)
        y = x.second[j]
        print(io,"$y \\[")
        Similitude.latexgroup(io,Similitude.evaldim(y))
        print(io,"\\]")
        j<length(x.second) && print(io,", ")
    end
    Values(str1,String(take!(io)))
end

function latexquantity(q::Group)
    io = IOBuffer()
    FieldAlgebra.special_print(io,product(q))
    print(io, " \\[\\mathbb{1}\\]")
    str1 = String(take!(io))
    FieldAlgebra.latexgroup_pre(io,q,FieldAlgebra.latext(q),'1')
    Values(str1,String(take!(io)),"Universe")
end
function latexquantity(q::LogGroup)
    io = IOBuffer()
    FieldAlgebra.special_print(io,product(q))
    print(io, " \\[\\mathbb{1}\\]")
    Values(String(take!(io)),"","Universe") # skip latexgroup_pre for log (sackur tetrode)
end
function latexquantity(q::Quantity{U}) where U
    v = 𝟏*q.v
    io = IOBuffer()
    FieldAlgebra.special_print(io,product(q.v))
    print(io, " \\[")
    latexgroup(io,normal(U)(dimensions(q)),U)
    print(io, "\\]")
    str1 = String(take!(io))
    FieldAlgebra.latexgroup_pre(io,q.v,FieldAlgebra.latext(v),'1')
    Values(str1,String(take!(io)),unitname(U))
end
function latexquantity(q::Similitude.ConvertUnit{U,S})  where {U,S}
    io = IOBuffer()
    rat = ratio(dimensions(q),U,S)
    FieldAlgebra.special_print(io,product(rat))
    d = convertdim(dimensions(q),U,S)
    print(io, " \\[")
    latexgroup(io,normal(S)(d),S)
    print(io, "\\]/\\[")
    latexgroup(io,normal(U)(d),U)
    print(io, "\\]")
    str1 = String(take!(io))
    FieldAlgebra.latexgroup_pre(io,rat,FieldAlgebra.latext(rat),'1')
    Values(str1,String(take!(io)),"$(unitname(U)) -> $(unitname(S))")
end

function latexdimensions(D,U)
    io = IOBuffer()
    latexgroup(io,U(D),U)
    String(take!(io))
end
latexdimensions(D,U::Tuple) = latexdimensions.(Ref(D),U)
latexdimensions(D::Vector,U::Tuple) = latexdimensions.(D,Ref(U))

if haskey(ENV,"UNITDOCS")
println("Similitude: documenting physics dimensions")
include("$dir/kinematicdocs.jl")
include("$dir/electromagneticdocs.jl")
include("$dir/thermodynamicdocs.jl")
println("Similitude: documenting physics constants")
include("$dir/physicsdocs.jl")
println("Similitude: documenting UnitSystems")
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
