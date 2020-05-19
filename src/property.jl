export SingleDofRSFProperty, RateStateQuasiDynamicProperty

## Property interface
import Base.fieldnames
import Base.==

abstract type AbstractProperty end

# https://github.com/jw3126/Setfield.jl
@doc raw"""
System property for single degree of freedom under rate-state friction.

## Fields
- `a`: contrib from velocity
- `b`: contrib from state
- `L`: critical distance
- `k`: spring stiffness
- `σ`: effective normal stress
- `η`: radiation damping
- `vpl`: plate rate
- `f0` = 0.6: ref. frictional coeff
- `v0` = 1e-6: ref. velocity
"""
@with_kw mutable struct SingleDofRSFProperty{T<:Real} <: AbstractProperty
    a::T # contrib from velocity
    b::T # contrib from state
    L::T # critical distance
    k::T # spring stiffness
    σ::T # effective normal stress
    η::T # radiation damping
    vpl::T # plate rate
    f0::T = 0.6 # ref. frictional coeff
    v0::T = 1e-6 # ref. velocity

    @assert a > 0
    @assert b > 0
    @assert L > 0
    @assert k > 0
    @assert σ > 0
    @assert η ≥ 0
    @assert vpl > 0
    @assert f0 > 0
    @assert v0 > 0
end

@doc raw"""
System property for multiple fault patches under rate-state friction.

## Fields
- `a`: contrib from velocity
- `b`: contrib from state
- `L`: critical distance
- `σ`: effective normal stress
- `η`: radiation damping
- `vpl`: plate rate
- `f0` = 0.6: ref. frictional coeff
- `v0` = 1e-6: ref. velocity
"""
@with_kw struct RateStateQuasiDynamicProperty{T<:Real, U<:AbstractVecOrMat} <: AbstractProperty
    a::U # contrib from velocity
    b::U # contrib from state
    L::U # critical distance
    σ::U # effective normal stress
    η::T # radiation damping
    vpl::T # plate rate
    f0::T = 0.6 # ref. frictional coeff
    v0::T = 1e-6 # ref. velocity

    @assert size(a) == size(b)
    @assert size(b) == size(L)
    @assert size(L) == size(σ)
    @assert f0 > 0
    @assert v0 > 0
    @assert η > 0
    @assert vpl > 0
end

const prop_field_names = Dict(
    :SingleDofRSFProperty => ("a", "b", "L", "k", "σ", "η", "vpl", "f0", "v0"),
    :RateStateQuasiDynamicProperty => ("a", "b", "L", "σ", "η", "vpl", "f0", "v0"),
    )

for (nn, fn) in prop_field_names
    @eval begin
        fieldnames(p::$(nn)) = $(fn)
        description(p::$(nn)) = String($(QuoteNode(nn)))
    end
end

function Base.:(==)(p1::P, p2::P) where P<:AbstractProperty
    reduce(&, [getfield(p1, Symbol(name)) == getfield(p2, Symbol(name)) for name in fieldnames(p1)])
end

## shortcut function
friction(flf::FrictionLawForm, v::T, θ::T, p::SingleDofRSFProperty) where T = friction(flf, v, θ, p.a, p.b, p.L, p.f0, p.v0)
friction(flf::FrictionLawForm, u::AbstractVecOrMat{T}, p::SingleDofRSFProperty) where T<:Real = friction(flf, u[1], u[2], p)
friction(flf::FrictionLawForm, u::AbstractArray{T}, p::SingleDofRSFProperty) where T<:AbstractVecOrMat = friction.(Ref(flf), u, Ref(p))
