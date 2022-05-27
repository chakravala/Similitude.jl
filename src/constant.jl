
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

include("$dir/constant.jl")

@pure factorize(x,a,b,c,d,e,f,g,h) = Similitude.factorize(x,param(a),param(b),param(c),param(d),param(e),param(f),param(g),param(h))

@pure Constant(N::Number) = 𝟏*N
@pure Constant(D::Dimension) = D
@pure Constant(::typeof(MathConstants.φ)) = φ
@pure Constant(::typeof(MathConstants.γ)) = Constant(valueat(35,vals))
@pure Constant(::typeof(ℯ)) = Constant(valueat(36,vals))
@pure Constant(::typeof(π)) = τ/𝟐
@pure Constant(::typeof(exp)) = Constant(ℯ)
@pure Constant(N::Float64) = Constant(factorize(N,τ,𝟐,𝟑,𝟓,𝟕,𝟏𝟏,𝟏𝟗,𝟒𝟑))
@pure Constant(N::Int) = Constant(factorize(N,τ,𝟐,𝟑,𝟓,𝟕,𝟏𝟏,𝟏𝟗,𝟒𝟑))
@pure Constant(N::AbelianGroup) = Constant{N}()

#printone(io::IO,::Val{vals}) = print(io, '𝟏')
Base.show(io::IO,x::Constant{N}) where N  = (showgroup(io,N,Natural,'𝟏'); print(io, " = ", constant(N)))

@pure constant(::Constant{N}) where N = constant(N)
constant(d::LogGroup{B},C=UnitSystems.Universe) where B = log(B,constant(value(d),C))
constant(d::ExpGroup{B},C=UnitSystems.Universe) where B = B^constant(value(d),C)

UnitSystems.unit(x::Constant,y=1) = x
@pure promoteint(v::Constant) = isone(v) ? 1 : v

Base.:^(::UnitSystems.Constant{A},::Constant{B}) where {A,B} = Constant{A^B}()
Base.:*(::UnitSystems.Constant{A},b::Constant) where A = Constant(A)*b
Base.:*(a::Constant,::UnitSystems.Constant{B}) where B = a*Constant(B)
Base.:/(::UnitSystems.Constant{A},b::Constant) where A = Constant(A)/b
Base.:/(a::Constant,::UnitSystems.Constant{B}) where B = a/Constant(B)

#Base.:*(a::Constant{N},b::Dimension{D}) where {N,D} = Dimension{N*D}()
#Base.:*(a::Dimension{D},b::Constant{N}) where {D,N} = Dimension{D*N}()
Base.:/(a::Constant,b::Dimension) = a*inv(b)
Base.:/(a::Dimension,b::Constant) = a*inv(b)

coefprod(a::Constant,b) = a*Constant(b)
coefprod(a,b::Constant) = Constant(a)*b
coefprod(a::Constant,b::Constant) = a*b

Base.:+(a::Constant,b::Similitude.Quantity{D,U}) where {D,U} = U(D)==𝟙 ? Similitude.Quantity{D,U}(a+b.v) : throw(error("$(U(D)) ≠ 𝟙 "))
Base.:+(a::Similitude.Quantity{D,U},b::Constant) where {D,U} = U(D)==𝟙 ? Similitude.Quantity{D,U}(a.v+b) : throw(error("$(U(D)) ≠ 𝟙 "))
Base.:-(a::Constant,b::Similitude.Quantity{D,U}) where {D,U} = U(D)==𝟙 ? Similitude.Quantity{D,U}(a-b.v) : throw(error("$(U(D)) ≠ 𝟙 "))
Base.:-(a::Similitude.Quantity{D,U},b::Constant) where {D,U} = U(D)==𝟙 ? Similitude.Quantity{D,U}(a.v-b) : throw(error("$(U(D)) ≠ 𝟙 "))
Base.:*(a::Constant,b::Similitude.Quantity{D,U}) where {D,U} = Similitude.Quantity{D,U}(a*b.v)
Base.:*(a::Similitude.Quantity{D,U},b::Constant) where {D,U} = Similitude.Quantity{D,U}(a.v*b)
Base.:/(a::Constant,b::Similitude.Quantity{D,U}) where {D,U} = Similitude.Quantity{inv(D),U}(a/b.v)
Base.:/(a::Similitude.Quantity{D,U},b::Constant) where {D,U} = Similitude.Quantity{D,U}(a.v/b)

phys(j,k=vals) = Constant(valueat(j,k))

for i ∈ 1:vals-11
    @eval begin
        export $(Symbol(basis[i]))
        const $(Symbol(basis[i])) = Constant(valueat($i,vals))
    end
end

export factorize

const golden = Constant(valueat(34,vals))
const eulergamma = Constant(valueat(35,vals))
const tau = Constant(valueat(37,vals))
const 𝟏 = Constant(valueat(0,vals))
const two = Constant(valueat(38,vals))
const three = Constant(valueat(39,vals))
const five = Constant(valueat(40,vals))
const seven = Constant(valueat(41,vals))
const eleven = Constant(valueat(42,vals))
const nineteen = Constant(valueat(43,vals))
const fourtythree = Constant(valueat(44,vals))
const 𝟐,𝟑,𝟓,𝟕,𝟏𝟏,𝟏𝟗,𝟒𝟑 = two,three,five,seven,eleven,nineteen,fourtythree
const zetta,zepto,yotta,yocto = (𝟐*𝟓)^21, (𝟐*𝟓)^-21, (𝟐*𝟓)^24, (𝟐*𝟓)^-24
const αinv,φ,τ = inv(α),golden,tau
const RK1990,KJ1990 = RK90,KJ90
const RK2014,KJ2014 = RK,KJ
