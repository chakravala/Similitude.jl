
#   This file is part of Similitude.jl
#   It is licensed under the AGPL license
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

function unitdim(D,U,S,L)
    addgroup(D,U,S)
    addlatex(D,U,L)
end

"""
    @unitdim(D,U,S,L)

Specify print `S::String` and LaTeX `L::String` for derived dimension `D` in `U::UnitSystem`.
```Julia
@unitdim magneticflux Gauss "Mx" "\\text{Mx}"
@unitdim magneticfluxdensity Gauss "G" "\\text{G}"
@unitdim magneticfield Gauss "Oe" "\\text{Oe}"
@unitdim frequency Metric "Hz" "\\text{Hz}"
@unitdim force Metric "N" "\\text{N}"
@unitdim pressure Metric "Pa" "\\text{Pa}"
@unitdim energy Metric "J" "\\text{J}"
@unitdim power Metric "W" "\\text{W}"
@unitdim mass British "slug" "\\text{slug}"
@unitdim force FPS "pdl" "\\text{pdl}"
```
These standard examples are some of the built-in defaults.
"""
macro unitdim(D, U, S, L)
    return :(unitdim($U($D),normal($U),$S,$L))
    #=quote
        Similitude.showgroup(io::IO,::typeof(Constant($U($D))),::typeof(normal($U))) = print(io,$S)
        Similitude.latexgroup(io::IO,::typeof(Constant($U($D))),::typeof(normal($U))) = print(io,$L)
    end=#
end

for U ∈ (:Engineering,:Gravitational)
    @eval begin
        @unitdim frequency $U "Hz" "\\text{Hz}"
        @unitdim frequencydrift $U "Hz⋅s⁻¹" "\\text{Hz} \\cdot \\text{s}^{-1}"
        @unitdim photonirradiance $U "Hz⋅m⁻²" "\\text{Hz} \\cdot \\text{m}^{-2}"
        @unitdim illuminance $U "lx" "\\text{lx}"
        @unitdim luminousexposure $U "lx⋅s" "\\text{lx} \\cdot \\text{s}"
        unitdim(luminance,normal($U),"nt","\\text{nt}")
        #showgroup(io::IO,::typeof(Constant(luminance)),::typeof(normal($U))) = print(io,"nt")
        #latexgroup(io::IO,::typeof(Constant(luminance)),::typeof(normal($U))) = print(io,"\\text{nt}")
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
        unitdim(luminance,normal($U),"nt","\\text{nt}")
        #Similitude.showgroup(io::IO,::typeof(Constant(luminance)),::typeof(normal($U))) = print(io,"nt")
        #Similitude.latexgroup(io::IO,::typeof(Constant(luminance)),::typeof(normal($U))) = print(io,"\\text{nt}")
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
                @unitdim radiantintensity $U $("W⋅$(u)⁻²") $("\\text{W} \\cdot $(uu)^{-2}")
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
unitdim(luminance,normal(Meridian),"ent","\\text{ent}")
#showgroup(io::IO,::typeof(Constant(luminance)),::typeof(normal(Meridian))) = print(io,"ent")
#latexgroup(io::IO,::typeof(Constant(luminance)),::typeof(normal(Meridian))) = print(io,"\\text{ent}")

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
        unitdim(luminance,normal($U),"sb","\\text{sb}")
        #showgroup(io::IO,::typeof(Constant(luminance)),::typeof(normal($U))) = print(io,"sb")
        #latexgroup(io::IO,::typeof(Constant(luminance)),::typeof(normal($U))) = print(io,"\\text{sb}")

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
@unitdim dimensions(molarmass) Gravitational "hyl⋅mol⁻¹" "\\text{hyl} \\cdot \\text{mol}^{-1}"
@unitdim dimensions(molarmass) British "slug⋅slug-mol⁻¹" "\\text{slug} \\cdot \\text{slug-mol}^{-1}"
@unitdim dimensions(molarmass) IPS "slinch-slinch-mol⁻¹" "\\text{slinch} \\cdot \\text{slinch-mol}^{-1}"
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
