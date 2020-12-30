using Quaycle
using GmshTools
using HDF5
using ProgressMeter

##
function events_index(t::T, maxv::T, thrd::U=1e-3; dura::U=1e3, shift::I=1, stride::I=1) where {T, U, I}
    lent = length(t)
    isevent = false
    idxb, idxe, idxe_unshift = Int[], Int[], Int[]
    for i ∈ eachindex(t)
        if maxv[i] ≥ thrd && !isevent
            isevent = true
            push!(idxb, (i-1)*stride+shift)
            _duraidx = findnearest(t, t[i] + dura)
            _idxe = findlast(x -> x ≥ thrd, view(maxv, i: _duraidx)) + i - 1
            push!(idxe_unshift, _idxe)
            push!(idxe, (_idxe-1)*stride+shift)
        elseif !isempty(idxe) && i > idxe_unshift[end]
            isevent = false
        end
    end
    return idxb, idxe
end

function findnearest(a::A, v::T; kwargs...) where {A, T}
    lena = length(a)
    i1 = searchsortedfirst(a, v; kwargs...)
    i1 > lena && return lena
    i1 == 1 && return 1
    abs(a[i1] - v) > abs(a[i1-1] - v) ? i1 - 1 : i1
end

function rupture_initial(f::AbstractString, idxb::T, patch::AbstractArray) where T<:Integer
    @assert ishdf5(f) "Not an HDF5 file."
    h5open(f, "r") do fid
        vd = d_open(fid, "v")
        vi = vd[axes(vd)[1:end-1]..., idxb]
        vi .*= patch
        argmax(vi)
    end
end

rupture_initial(f::AbstractString, idxb::AbstractVector, patch::AbstractArray) = map(x -> rupture_initial(f, x, patch), idxb)

function accumulated_slip(δd::HDF5Dataset, idxb::T, idxe::T) where T<:Integer
    δax = axes(δd)[1:end-1]
    return δd[δax...,idxe] - δd[δax...,idxb]
end

accumulated_slip(f::S, idxb::T, idxe::T; δstr::S="δ", kwargs...) where {T<:Integer, S<:AbstractString} =
    h5open(f, "r") do fid
        accumulated_slip(d_open(fid, δstr), idxb, idxe; kwargs...)
    end

accumulated_slip(f::AbstractString, idxb::T, idxe::T; kwargs...) where T<:AbstractVector =
    map((b, e) -> accumulated_slip(f, b, e; kwargs...), idxb, idxe)

slip_magnitude(s::T, ds::T, μ::T=3e10) where T = (log10(μ * s * ds) - 9.1) / 1.5
slip_magnitude(s::AbstractVector, ds::T, μ::T=3e10) where T = map(x -> slip_magnitude(x, ds, μ), s)

##
function scan_output_01(filename; _stride=1, basedir="", edge=1e3)
    fname = joinpath(dirname(@__DIR__), "out", basedir, filename)
    include(joinpath(@__DIR__, "s01-domain.jl"))
    @info "Reading: " * fname

    f = h5open(fname, "r")
    td = d_open(f, "t")
    vd = d_open(f, "v")
    lent = length(td)
    iterindex = 1: _stride: lent # no shift
    ti = td[iterindex]
    @info "Last year: " td[lent]/365/86400
    patch_A = (@. -40e3 ≤ mf.x ≤ -edge) .* (@. -7e3 ≤ mf.z ≤ 0e3)'
    patch_B = (@. edge ≤ mf.x ≤ 40e3) .* (@. -7e3 ≤ mf.z ≤ 0e3)'
    patch_F = trues(mf.nx, mf.nξ)
    patches = [patch_A, patch_B, patch_F]
    maxvs = [Vector{eltype(td[1])}(undef, length(iterindex)) for _ in 1: 3]
    vsz = size(vd)
    tmps = [Array{eltype(td[1])}(undef, vsz[1:end-1]...) for _ in 1: 3]
    @info "Time steps: " * string(lent)
    iterindex_ = collect(iterindex)
    leniter = length(iterindex_)
    @inbounds for i ∈ eachindex(iterindex_)
        i % 100 == 0 && @info i/leniter
    	vi = vd[:, :, iterindex_[i]]
    	for j ∈ 1: 3
            @. tmps[j] = vi * patches[j]
            maxvs[j][i] = maximum(tmps[j])
    	end
    end

    ## 2e5 seconds for upper limit of seismic duration
    ixb1, ixe1 = events_index(ti, maxvs[1], 1e-3; dura=2e5, stride=_stride)
    ixb2, ixe2 = events_index(ti, maxvs[2], 1e-3; dura=2e5, stride=_stride)

    ss1 = accumulated_slip(fname, ixb1, ixe1)
    ss2 = accumulated_slip(fname, ixb2, ixe2)
    sss1 = [sum(x .* patch_A) for x in ss1]
    sss2 = [sum(x .* patch_B) for x in ss2]
    mw1 = slip_magnitude(sss1, mf.Δx * mf.Δξ)
    mw2 = slip_magnitude(sss2, mf.Δx * mf.Δξ)
    ri1 = rupture_initial(fname, ixb1, patch_A)
    ri2 = rupture_initial(fname, ixb2, patch_B)
    riax = [mf.x[q.I[1]] for q in ri1]
    riaz = [mf.ξ[q.I[2]] for q in ri1]
    ribx = [mf.x[q.I[1]] for q in ri2]
    ribz = [mf.ξ[q.I[2]] for q in ri2]

    if isempty(ixb1) || isempty(ixe1)
        riax, riaz, mw1 = Float64[], Float64[], Float64[]
    elseif isempty(ixb2) || isempty(ixe2)
        ribx, ribz, mw2 = Float64[], Float64[], Float64[]
    end

    results = joinpath(dirname(@__DIR__), "out", basedir, "res" * basename(fname)[3:end])
    isfile(results) && rm(results)
    h5write(results, "maxva", maxvs[1])
    h5write(results, "maxvb", maxvs[2])
    h5write(results, "maxvf", maxvs[3])
    h5write(results, "t", td[:])
    h5write(results, "ti", ti)
    h5write(results, "riax", riax)
    h5write(results, "ribx", ribx)
    h5write(results, "riaz", riaz)
    h5write(results, "ribz", ribz)
    h5write(results, "mwa", mw1)
    h5write(results, "mwb", mw2)
    h5write(results, "ixba", ixb1)
    h5write(results, "ixbb", ixb2)
    h5write(results, "ixea", ixe1)
    h5write(results, "ixeb", ixe2)

    return results
end

function coulomb_stress(gf, δ)
    nx, ny = size(gf)[1:2]
    δa = cat(δ, zeros(nx-1, ny); dims=1)
    τka = zeros(ComplexF64, nx, ny)
    δka = rfft(δa, 1)
    for j ∈ 1: ny, l ∈ 1: ny, i ∈ 1: nx
        τka[i,j] += gf[i,j,l] * δka[i,l]
    end
    τa = irfft(τka, 2*nx-1, 1) # works for even number of `nx`
    τ = τa[1:nx,:]
end

function slip_ratio(fname, cutoff=:all; basedir="")
    f = joinpath(dirname(@__DIR__), "out", basedir, fname)
    include(joinpath(@__DIR__, "s01-domain.jl"))
    @info "Reading $(f)"

    t = h5read(f, "t")
    final_slip = h5read(f, "δ", (:,:,length(t)))
    if cutoff == :all
        cutoff = size(final_slip, 2)
    else
        cutoff = findnearest(mf.ξ, cutoff, rev=true)
    end
    @info "Cutoff $(mf.ξ[cutoff])"

    maxv = zeros(size(final_slip, 1), length(t)) # along strike
    for i in eachindex(t)
        i % 1000 == 0 && @info "Scan Maxv: $(i/length(t))"
        maxv[:, i] .= maximum(h5read(f, "v", (:, :, i)); dims=2) |> vec
    end

    seis_slip = zeros(size(final_slip))
    after_slip = zeros(size(final_slip))

    for i in 1: size(maxv, 1)
        @info "Compute sr index: $(i)"
        ixb_c, ixe_c = events_index(t, maxv[i, :], 1e-3; dura=2e5)
        ixb_a, ixe_a = events_index(t, maxv[i, :], 1e-6; dura=3e7)
        ixe_aa = map(x -> ixe_a[findfirst(>=(x), ixe_a)], ixe_c)
        @assert all(ixe_aa .≥ ixe_c)
        for (b, e) in zip(ixb_c, ixe_c)
            seis_slip[i, :] .+= h5read(f, "δ", (i,:,e)) - h5read(f, "δ", (i,:,b))
        end
        for (b, e) in zip(ixe_c, ixe_aa)
            after_slip[i, :] .+= h5read(f, "δ", (i,:,e)) - h5read(f, "δ", (i,:,b))
        end
    end

    seis_ratio = sum(seis_slip[:, 1: cutoff]; dims=2) ./ sum(final_slip[:, 1: cutoff]; dims=2) |> vec
    after_ratio = sum(after_slip[:, 1: cutoff]; dims=2) ./ sum(final_slip[:, 1: cutoff]; dims=2) |> vec

    results = joinpath(dirname(@__DIR__), "out", basedir, "sr" * basename(fname)[3:end])
    isfile(results) && rm(results)
    h5write(results, "sr", seis_ratio)
    h5write(results, "ar", after_ratio)
end
