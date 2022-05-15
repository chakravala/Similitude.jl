
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

include("$dir/constant.jl")

@pure factorize(x,a,b,c,d,e,f,g,h) = Similitude.factorize(x,param(a),param(b),param(c),param(d),param(e),param(f),param(g),param(h))

@pure Constant(D::Number) = 𝟏*D
@pure Constant(D::Float64) = Constant(factorize(D,τ,𝟐,𝟑,𝟓,𝟕,𝟏𝟏,𝟏𝟗,𝟒𝟑))
@pure Constant(D::Int) = Constant(factorize(D,τ,𝟐,𝟑,𝟓,𝟕,𝟏𝟏,𝟏𝟗,𝟒𝟑))
@pure Constant(D::Constant) = D
@pure Constant(D::AbelianGroup) = Constant{D}()

#printone(io::IO,::Val{vals}) = print(io, '𝟏')
Base.show(io::IO,x::Constant{D}) where D  = (showgroup(io,D,Natural,'𝟏'); print(io, " = ", constant(D)))

@pure constant(::Constant{D}) where D = constant(D)
constant(d::LogGroup{B},C=UnitSystems.Universe) where B = log(B,constant(value(d),C))
constant(d::ExpGroup{B},C=UnitSystems.Universe) where B = B^constant(value(d),C)

UnitSystems.unit(x::Constant,y=1) = x

Base.:*(a::UnitSystems.Constant{A},b::Constant{B}) where {A,B} = Constant{A*B}()
Base.:*(a::Constant{A},b::UnitSystems.Constant{B}) where {A,B} = Constant{A*B}()
Base.:/(a::UnitSystems.Constant{A},b::Constant{B}) where {A,B} = Constant{A/B}()
Base.:/(a::Constant{A},b::UnitSystems.Constant{B}) where {A,B} = Constant{A/B}()

Base.:+(a::Constant,b::Similitude.Quantity{D,U}) where {D,U} = U(D)==𝟙 ? Similitude.Quantity{D,U}(a+b.v) : throw(error("$(U(D)) ≠ 𝟙 "))
Base.:+(a::Similitude.Quantity{D,U},b::Constant) where {D,U} = U(D)==𝟙 ? Similitude.Quantity{D,U}(a.v+b) : throw(error("$(U(D)) ≠ 𝟙 "))
Base.:-(a::Constant,b::Similitude.Quantity{D,U}) where {D,U} = U(D)==𝟙 ? Similitude.Quantity{D,U}(a-b.v) : throw(error("$(U(D)) ≠ 𝟙 "))
Base.:-(a::Similitude.Quantity{D,U},b::Constant) where {D,U} = U(D)==𝟙 ? Similitude.Quantity{D,U}(a.v-b) : throw(error("$(U(D)) ≠ 𝟙 "))
Base.:*(a::Constant,b::Similitude.Quantity{D,U}) where {D,U} = Similitude.Quantity{D,U}(a*b.v)
Base.:*(a::Similitude.Quantity{D,U},b::Constant) where {D,U} = Similitude.Quantity{D,U}(a.v*b)
Base.:/(a::Constant,b::Similitude.Quantity{D,U}) where {D,U} = Similitude.Quantity{inv(D),U}(a/b.v)
Base.:/(a::Similitude.Quantity{D,U},b::Constant) where {D,U} = Similitude.Quantity{D,U}(a.v/b)

phys(j,k=vals) = Constant(valueat(j,k))

for i ∈ 1:vals-10
    @eval begin
        export $(Symbol(basis[i]))
        const $(Symbol(basis[i])) = Constant(valueat($i,vals))
    end
end

export factorize

const τ = Constant(valueat(37,vals))
const 𝟏 = Constant(valueat(0,vals))
const 𝟐 = Constant(valueat(38,vals))
const 𝟑 = Constant(valueat(39,vals))
const 𝟓 = Constant(valueat(40,vals))
const 𝟕 = Constant(valueat(41,vals))
const 𝟏𝟏 = Constant(valueat(42,vals))
const 𝟏𝟗 = Constant(valueat(43,vals))
const 𝟒𝟑 = Constant(valueat(44,vals))
const 𝟏𝟎 = 𝟐*𝟓
const zetta,zepto,yotta,yocto = (𝟐*𝟓)^21, (𝟐*𝟓)^-21, (𝟐*𝟓)^24, (𝟐*𝟓)^-24
const αinv = inv(α)
const RK1990,KJ1990 = RK90,KJ90
const RK2014,KJ2014 = RK,KJ
