
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

import UnitSystems: isconstant

using FieldConstants
import FieldConstants: Constant

@pure constant(N::Number) = 𝟏*N
@pure constant(D::Constant) = D
@pure constant(::typeof(MathConstants.φ)) = φ
@pure constant(::typeof(MathConstants.γ)) = Constant(valueat(35,vals,:Constants))
@pure constant(::typeof(ℯ)) = Constant(valueat(36,vals,:Constants))
@pure constant(::typeof(π)) = τ/𝟐
@pure constant(::typeof(exp)) = constant(ℯ)
@pure constant(N::Float64) = Constant(factorize(N,Val(:Constants)))
@pure constant(N::Int) = Constant(factorize(N,τ,Val(:Constants)))
@pure constant(N::AbelianGroup) = Constant{N}()

#constant(d::LogGroup{B},C=UnitSystems.Universe) where B = log(B,constant(value(d),C))
#constant(d::ExpGroup{B},C=UnitSystems.Universe) where B = B^constant(value(d),C)

export factorize

const golden = phys(34)
const eulergamma = phys(35)
const tau = phys(37)
const 𝟏 = phys(0)
const two = phys(38)
const three = phys(39)
const five = phys(40)
const seven = phys(41)
const eleven = phys(42)
const nineteen = phys(43)
const fourtythree = phys(44)
const 𝟐,𝟑,𝟓,𝟕,𝟏𝟏,𝟏𝟗,𝟒𝟑 = two,three,five,seven,eleven,nineteen,fourtythree
const zetta,zepto,yotta,yocto = (𝟐*𝟓)^21, (𝟐*𝟓)^-21, (𝟐*𝟓)^24, (𝟐*𝟓)^-24
const αinv = inv(α)
const RK1990,KJ1990 = RK90,KJ90
const RK2014,KJ2014 = RK,KJ
