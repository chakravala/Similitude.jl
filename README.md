# Similitude.jl

*Dimensions and Quantities for UnitSystems*

[![DOI](https://zenodo.org/badge/320717758.svg)](https://zenodo.org/badge/latestdoi/320717758)
[![Build status](https://ci.appveyor.com/api/projects/status/3fd8dauls91okw8q?svg=true)](https://ci.appveyor.com/project/chakravala/similitude-jl)

In the base `UnitSystems` package, simply `Float64` numbers are used to generate the group of `UnitSystem` constants.
However, in the `Similitude` package, instead constant numbers are used to generate an abstract multiplicative `Group`, which is only converted to a `Float64` value when appropriate or at compile time.

In `UnitSystems` and `Similitude`, the spectrum of `Group{:Constants}` values is generated by a group of 11 mathematical constants (7 Integer primes and 4 Irrational numbers) with 33 physical measurement definitions.
These are `𝟐`, `𝟑`, `𝟓`, `𝟕`, `𝟏𝟏`, `𝟏𝟗`, `𝟒𝟑`, `φ`, `γ`, `ℯ`, `τ`, `kB`, `NA`, `𝘩`, `𝘤`, `𝘦`, `Kcd`, `ΔνCs`, `R∞`, `α`, `μₑᵤ`, `μₚᵤ`, `ΩΛ`, `H0`, `g₀`, `aⱼ`, `au`, `ft`, `ftUS`, `lb`, `T₀`, `atm`, `inHg`, `RK90`, `KJ90`, `RK`, `KJ`, `Rᵤ2014`, `Ωᵢₜ`, `Vᵢₜ`, `kG`, `mP`, `GME`, `GMJ`.

Furthermore, in `Similitude` there is a dimension type which encodes the dimensional `Group{:USQ}` for the `Quantity` type using the same implementation principles as the `Group{:Constants}`.
This enables the unified usage of `Group` homomorphisms to transform `Quantity` algebra elements with varying numbers of dimensionless constants.

Originally, the Newtonian group used for `UnitSystems` would be made up of `force`, `mass`, `length`, `time` (or `F`, `M`, `L`, `T`).
Although `force` is typically thought of as a derived dimension when the reference `gravity` is taken to be dimensionless, `force` is actually considered a base dimension in general engineering `UnitSystem` foundations.
With the development of electricity and magnetism came an interest for an additional dimension called `charge` or `Q`.
When the thermodynamics of `entropy` became further developed, the `temperature` or `Θ` was introduced as another dimension.
In the field of chemistry, it became desirable to introduce another dimension of `molaramount` or `N` as fundamental.
To complete the existing International System of Quantities (ISQ) it is also necessary to consider `luminousflux` or `J` as a visual perception related dimension.
In order to resolve ambiguity with `solidangle` unit conversion, `angle` or `A` is explicitly tracked in the underlying dimension and `Group`.
However, this is yet insufficient to fully specify all the historical variations of `UnitSystem`, including the `EMU`, `ESU`, `Gauss` and `LorentzHeavise` specifications.
Therefore, there is also a dimension basis for `rationalization` (denoted `R`) and `lorentz` (denoted by `C⁻¹`).

In combination, all these required base dimension definitions are necessary in order to coherently implement unit conversion for `Quantity` elements.
Since the existing International System of Quantities (ISQ) is an insufficient definition for dimension, a new Unified System of Quantities (`USQ`) is being proposed here as composed of `force`, `mass`, `length`, `time`, `charge`, `temperature`, `molaramount`, `luminousflux`, `angle`, `rationalization`, and a `nonstandard` dimension (denoted by `F`, `M`, `L`, `T`, `Q`, `Θ`, `N`, `J`, `A`, `R`, `C`).

In aggregate, the `UnitSystem` data generated here constitutes a new universal standardization for dimensional analysis, which generalizes upon previous historical systems up to the 2019 redefinition and unifies them in a common `Universe`.
This enables a more precise and generalized standardization than the 2019 redefinition, which was comparatively limited in scope.
Specified default `UnitSystem` values are to be taken as a newly defined mutually-compatible recommended standard, verified to be consistent and coherent.

> In fact there is nothing transcendental about dimensions; the ultimate principle is precisely expressible (in Newton's terminology) as one of *similitude*, exact or approximate, to be tested by the rule that mere change in the magnitudes of the ordered scheme of units of measurement that is employed must not affect sensibly the forms of the equations that are the adequate expression of the underlying relations of the problem. (J.L., 1914)

Specifications for dimensional units are in the [UnitSystems.jl](https://github.com/chakravala/UnitSystems.jl) and [Similitude.jl](https://github.com/chakravala/Similitude.jl) and [MeasureSystems.jl](https://github.com/chakravala/MeasureSystems.jl) repositories.
The three packages are designed so that they can be interchanged with compatibility.
On its own `UnitSystems` is the fastest package, while `Similitude` (provides `Quantity` type) and `MeasureSystems` (introduces [Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl) uncertainty) build additional features on top of `UnitSystems` base defintions.
Defaults are shared: `Metric`, `SI2019`, `CODATA`, `Conventional`, `International`, `InternationalMean`, `MetricTurn`, `MetricGradian`, `MetricDegree`, `MetricArcminute`, `MetricArcsecond`, `Engineering`, `Gravitational`, `FPS`, `IPS`, `British`, `English`, `Survey`, `Gauss`, `LorentzHeaviside`, `EMU`, `ESU`, `IAU`, `IAUE`, `IAUJ`, `Hubble`, `Cosmological`, `CosmologicalQuantum`, `Meridian`, `Nautical`, `MPH`, `KKH`, `MTS`, `FFF`, `Planck`, `PlanckGauss`, `Stoney`, `Hartree`, `Rydberg`, `Schrodinger`, `Electronic`, `Natural`, `NaturalGauss`, `QCD`, `QCDGauss`, `QCDoriginal`.

```Julia
julia> using Similitude # or UnitSystems or MeasureSystems
```

An optional environment variable `ENV["SIMILITUDE"]` induces `UnitSystems.similitude()` to return `true`, giving flexibility for building dependencies whenever it is desirable to toggle usage between `UnitSystems` (default) and `Similitude` (requires environment variable specification). For example, in `MeasureSystems` and `Geophysics` this option is used to increase flexibility with variety in local compilation workflow.

A `UnitSystem` is a consistent set of dimensional values selected to accomodate a particular use case or standardization.
It is possible to convert derived physical quantities from any `UnitSystem` specification into any other using accurate values.
Eleven fundamental constants `kB`, `ħ`, `𝘤`, `μ₀`, `mₑ`, `Mᵤ`, `Kcd`, `θ`, `λ`, `αL`, `g₀` are used to govern a specific unit system consistent scaling.
These are the constants `boltzmann`, `planckreduced`, `lightspeed`, `vacuumpermeability`, `electronmass`, `molarmass`, `luminousefficacy`, `angle`, `rationalization`, `lorentz`, and `gravity`.
Different choices of natural units or physical measurements result in a variety of unit systems for many purposes.

Main documentation is at https://geophysics.crucialflow.com/dev/unitsystems

Historically, older electromagnetic unit systems also relied on a `rationalization` constant `λ` and a `lorentz` force proportionality constant `αL`.
In most unit systems these extra constants have a value of `1` unless specified.

```Julia
    UnitSystem{kB, ħ, 𝘤, μ₀, mₑ, Mᵤ, (Kcd, θ, λ, αL, g₀, ...)}
```

Fundamental constants of physics are: `kB` Boltzmann's constant, `ħ` reduced Planck's constant, `𝘤` speed of light, `μ₀` vacuum permeability, `mₑ` electron rest mass, `Mᵤ` molar mass, `Kcd` luminous efficacy, `θ` angle measure, `λ` Gauss rationalization, `αL` Lorentz's constant, and `g₀` gravitational force reference.
Primarily the `Metric` SI unit system is used in addition to the historic `English` engineering unit system.
These constants induce derived values for `avogadro`, `boltzmann`, `molargas`, `planck`, `planckreduced`, `lightspeed`, `planckmass`, `dalton`, `protonmass`, `electronmass`, `newton`, `einstein`, `vacuumpermeability`, `vacuumpermittivity`, `electrostatic`, and
additional constants `molarmass`, `luminousefficacy`, `gravity`, `angle`, `turn`, `spat`, `stefan`, `radiationdensity`, `magnetostatic`, `lorentz`, `biotsavart`, `rationalization`, `vacuumimpedance`, `elementarycharge`, `magneton`, `conductancequantum`, `faraday`, `magneticfluxquantum`, `josephson`, `klitzing`, `hartree`, `rydberg`, `bohr`.

Physics constant documentation is at https://geophysics.crucialflow.com/dev/constants

Standardized unit/derived quantities are `hyperfine`, `loschmidt`, `wienwavelength`, `wienfrequency`, `mechanicalheat`, `eddington`, `solarmass`, `jupitermass`, `earthmass`, `lunarmass`, `earthradius`, `greatcircle`, `radarmile`, `hubble`, `cosmological`, `radian`, `steradian`, `degree`, `squaredegree`, `gradian`, `arcminute`, `arcsecond`, `second`, `minute`, `hour`, `day`, `gaussianmonth`, `siderealmonth`, `synodicmonth`, `year`, `gaussianyear`, `siderealyear`, `jovianyear`, `angstrom`, `inch`, `foot`, `surveyfoot`, `yard`, `meter`, `earthmeter`, `mile`, `statutemile`, `meridianmile`, `admiraltymile`, `nauticalmile`, `lunardistance`, `astronomicalunit`, `jupiterdistance`, `lightyear`, `parsec`, `bubnoff`, `ips`, `fps`, `fpm`, `ms`, `kmh`, `mph`, `knot`, `mps`, `barn`, `hectare`, `acre`, `surveyacre`, `liter`, `gallon`, `quart`, `pint`, `cup`, `fluidounce`, `teaspoon`, `tablespoon`, `grain`, `gram`, `earthgram`, `kilogram`, `tonne`, `ton`, `pound`, `ounce`, `slug`, `slinch`, `hyl`, `dyne`, `newton`, `poundal`, `poundforce`, `kilopond`, `psi`, `pascal`, `bar`, `barye`, `technicalatmosphere`, `atmosphere`, `inchmercury`, `torr`, `electronvolt`, `erg`, `joule`, `footpound`, `calorie`, `kilocalorie`, `meancalorie`, `earthcalorie`, `thermalunit`, `gasgallon`, `tontnt`, `watt`, `horsepower`, `horsepowerwatt`, `horsepowermetric`, `electricalhorsepower`, `tonsrefrigeration`, `boilerhorsepower`, `coulomb`, `earthcoulomb`, `ampere`, `volt`, `henry`, `ohm`, `siemens`, `farad`, `weber`, `tesla`, `abcoulomb`, `abampere`, `abvolt`, `abhenry`, `abohm`, `abmho`, `abfarad`, `maxwell`, `gauss`, `oersted`, `gilbert`, `statcoulomb`, `statampere`, `statvolt`, `stathenry`, `statohm`, `statmho`, `statfarad`, `statweber`, `stattesla`, `kelvin`, `rankine`, `celsius`, `fahrenheit`, `sealevel`, `boiling`, `mole`, `earthmole`, `poundmole`, `slugmole`, `slinchmole`, `katal`, `amagat`, `lumen`, `candela`, `lux`, `phot`, `footcandle`, `nit`, `apostilb`, `stilb`, `lambert`, `footlambert`, `bril`, `neper`, `bel`, `decibel`, `hertz`, `apm`, `rpm`, `kayser`, `diopter`, `gforce`, `galileo`, `eotvos`, `darcy`, `poise`, `reyn`, `stokes`, `rayl`, `mpge`, `langley`, `jansky`, `solarflux`, `curie`, `gray`, `roentgen`, `rem`.

Standard physics units are at https://geophysics.crucialflow.com/dev/units

Additional reference `UnitSystem` variants: `EMU`, `ESU`, `Gauss`, `LorentzHeaviside`, `SI2019`, `SI1976`, `CODATA`, `Conventional`, `International`, `InternationalMean`, `Engineering`, `Gravitational`, `IAU`, `IAUE`, `IAUJ`, `FPS`, `IPS`, `British`, `Survey`, `Hubble`, `Cosmological`, `CosmologicalQuantum`, `Meridian`, `Nautical`, `MPH`, `KKH`, `MTS`, `FFF`; and natural atomic units based on gravitational `coupling` and `finestructure` constant (`Planck`, `PlanckGauss`, `Stoney`, `Hartree`, `Rydberg`, `Schrodinger`, `Electronic`, `Natural`, `NaturalGauss`, `QCD`, `QCDGauss`, and `QCDoriginal`).

Unit conversion documentation is at https://geophysics.crucialflow.com/dev/convert

**Derived Unit conversions:**

Mechanics: `angle`, `solidangle`, `time`, `angulartime`, `length`, `angularlength`, `area`, `angulararea`, `volume`, `wavenumber`, `angularwavenumber`, `fuelefficiency`, `numberdensity`, `frequency`, `angularfrequency`, `frequencydrift`, `stagnance`, `speed`, `acceleration`, `jerk`, `snap`, `crackle`, `pop`, `volumeflow`, `etendue`, `photonintensity`, `photonirradiance`, `photonradiance`,
`inertia`, `mass`, `massflow`, `lineardensity`, `areadensity`, `density`, `specificweight`, `specificvolume`, `force`, `specificforce`, `gravityforce`, `pressure`, `compressibility`, `viscosity`, `diffusivity`, `rotationalinertia`, `impulse`, `momentum`, `angularmomentum`, `yank`, `energy`, `specificenergy`, `action`, `fluence`, `power`, `powerdensity`, `irradiance`, `radiance`, `radiantintensity`, `spectralflux`, `spectralexposure`, `soundexposure`, `impedance`, `specificimpedance`, `admittance`, `compliance`, `inertance`;
Electromagnetics: `charge`, `chargedensity`, `linearchargedensity`, `exposure`, `mobility`, `current`, `currentdensity`, `resistance`, `conductance`, `resistivity`, `conductivity`, `capacitance`, `inductance`, `reluctance`, `permeance`, `permittivity`, `permeability`, `susceptibility`, `specificsusceptibility`, `demagnetizingfactor`, `vectorpotential`, `electricpotential`, `magneticpotential`, `electricfield`, `magneticfield`, `electricflux`, `magneticflux`, `electricdisplacement`, `magneticfluxdensity`, `electricdipolemoment`, `magneticdipolemoment`, `electricpolarizability`, `magneticpolarizability`, `magneticmoment`, `specificmagnetization`, `polestrength`;
Thermodynamics: `temperature`, `entropy`, `specificentropy`, `volumeheatcapacity`, `thermalconductivity`, `thermalconductance`, `thermalresistivity`, `thermalresistance`, `thermalexpansion`, `lapserate`,
`molarmass`, `molality`, `mole`, `molarity`, `molarvolume`, `molarentropy`, `molarenergy`, `molarconductivity`, `molarsusceptibility`, `catalysis`, `specificity`, `diffusionflux`,
`luminousflux`, `luminousintensity`, `luminance`, `illuminance`, `luminousenergy`, `luminousexposure`, `luminousefficacy`.

Other similar packages include [UnitSystems.jl](https://github.com/chakravala/UnitSystems.jl), [MeasureSystems.jl](https://github.com/chakravala/MeasureSystems.jl), [PhysicalConstants.jl](https://github.com/JuliaPhysics/PhysicalConstants.jl), [MathPhysicalConstants.jl](https://github.com/LaGuer/MathPhysicalConstants.jl), [Unitful.jl](https://github.com/PainterQubits/Unitful.jl.git), [UnitfulUS.jl](https://github.com/PainterQubits/UnitfulUS.jl), [UnitfulAstro.jl](https://github.com/JuliaAstro/UnitfulAstro.jl), [UnitfulAtomic.jl](https://github.com/sostock/UnitfulAtomic.jl), [NaturallyUnitful.jl](https://github.com/MasonProtter/NaturallyUnitful.jl), and [UnitfulMoles.jl](https://github.com/rafaqz/UnitfulMoles.jl).
