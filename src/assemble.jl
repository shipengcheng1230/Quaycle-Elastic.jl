export assemble

"""
    assemble(p::SingleDofRSFProperty, u0::AbstractArray, tspan::NTuple;
        flf::FrictionLawForm=CForm(), se::StateEvolutionLaw=DieterichStateLaw())

Assemble the `ODEProblem` for single degree of freedom system.

## Arguments
- `p::SingleDofRSFProperty`: all system properties
- `u0::AbstractArray`: initial condition
- `tspan::NTuple`: time span for simulation
- `flf::FrictionLawForm`: form of friction law, either [`CForm`](@ref) or [`RForm`](@ref)
- `se::StateEvolutionLaw`: state evolutional law, see [`StateEvolutionLaw`](@ref)
"""
function assemble(p::SingleDofRSFProperty, u0::AbstractArray, tspan::NTuple;
    flf::FrictionLawForm=CForm(), se::StateEvolutionLaw=DieterichStateLaw()) where T<:Real
    (typeof(flf) == RForm && p.η ≈ 0) && @warn "Regularized form requires nonzero `η` to avoid `Inf` in dv/dt."
    op! = (du, u, p, t) -> du .= dv_dθ_dt(flf, se, u[1], u[2], p.a, p.b, p.L, p.k, p.σ, p.η, p.vpl, p.f0, p.v0)
    return ODEProblem(op!, u0, tspan, p)
end

"""
    assemble(gf::AbstractArray, p::RateStateQuasiDynamicProperty,
        u0::AbstractArray, tspan::NTuple{2};
        flf::FrictionLawForm=RForm(), se::StateEvolutionLaw=DieterichStateLaw(), kwargs...)

Assemble the `ODEProblem` for elastic fault.

##  Extra Arguments
- `gf::AbstractArray`: green's function associated with `fs.mesh` and `p.λ` & `p.μ`
- `p::RateStateQuasiDynamicProperty`: all system properties
- `u0::ArrayPartition`: initial condition. By rule of order in this package:
    1. *velocity*
    2. *state*


- `tspan::NTuple`: time span for simulation
- `flf::FrictionLawForm`: form of friction law, either [`CForm`](@ref) or [`RForm`](@ref)
- `se::StateEvolutionLaw`: state evolutional law, see [`StateEvolutionLaw`](@ref)
"""
function assemble(
    gf::AbstractArray, p::RateStateQuasiDynamicProperty, u0::ArrayPartition, tspan::NTuple{2};
    flf::FrictionLawForm=RForm(), se::StateEvolutionLaw=DieterichStateLaw(),
    )
    _assemble(gf, p, u0, tspan, flf, se)
end

@inline function _assemble(gf, p::AbstractProperty, u0, tspan, flf, se)
    alloc = gen_alloc(gf)::AbstractAllocation
    f! = (du, u, p, t) -> ∂u∂t(du, u, p, alloc, gf, flf, se)
    return ODEProblem(f!, u0, tspan, p)
end
