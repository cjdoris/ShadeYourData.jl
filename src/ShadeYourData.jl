module ShadeYourData

import Base.Threads: @threads
import Makie: Makie, Heatmap, Image, Axis, Theme, FRect3D, heatmap!, image!, lift, @recipe, default_theme, (..), observe_changes

struct Canvas{T}
    xmin :: T
    xmax :: T
    xsize :: Int
    ymin :: T
    ymax :: T
    ysize :: Int
end

Base.size(c::Canvas) = (c.xsize, c.ysize)
xlims(c::Canvas) = (c.xmin, c.xmax)
ylims(c::Canvas) = (c.ymin, c.ymax)
Base.:(==)(a::Canvas, b::Canvas) = size(a)==size(b) && xlims(a)==xlims(b) && ylims(a)==ylims(b)

abstract type AggOp end

update(a::AggOp, x, args...) = merge(a, x, embed(a, args...))

struct AggCount{T} <: AggOp end
AggCount() = AggCount{Int}()
null(::AggCount{T}) where {T} = zero(T)
embed(::AggCount{T}) where {T} = oneunit(T)
merge(::AggCount{T}, x::T, y::T) where {T} = x + y
value(::AggCount{T}, x::T) where {T} = x

struct AggAny <: AggOp end
null(::AggAny) = false
embed(::AggAny) = true
merge(::AggAny, x::Bool, y::Bool) = x | y
value(::AggAny, x::Bool) = x

struct AggSum{T} <: AggOp end
AggSum() = AggSum{Float64}()
null(::AggSum{T}) where {T} = zero(T)
embed(::AggSum{T}, x) where {T} = convert(T, x)
merge(::AggSum{T}, x::T, y::T) where {T} = x + y
value(::AggSum{T}, x::T) where {T} = x

struct AggMean{T} <: AggOp end
AggMean() = AggMean{Float64}()
null(::AggMean{T}) where {T} = (zero(T), zero(T))
embed(::AggMean{T}, x) where {T} = (convert(T,x), oneunit(T))
merge(::AggMean{T}, x::Tuple{T,T}, y::Tuple{T,T}) where {T} = (x[1]+y[1], x[2]+y[2])
value(::AggMean{T}, x::Tuple{T,T}) where {T} = float(x[1]) / float(x[2])

abstract type AggMethod end

struct AggSerial <: AggMethod end
struct AggThreads <: AggMethod end

aggregate(c::Canvas, xs, ys, zs...; op::AggOp=AggCount(), method::AggMethod=AggSerial()) =
    _aggregate(c, op, method, xs, ys, zs...)

function _aggregate(c::Canvas, op::AggOp, method::AggSerial, xs, ys, zs...)
    xmin, xmax = xlims(c)
    ymin, ymax = ylims(c)
    xsize, ysize = size(c)
    xmax > xmin || error("require xmax > xmin")
    ymax > ymin || error("require ymax > ymin")
    xscale = xsize / (xmax - xmin)
    yscale = ysize / (ymax - ymin)
    out = fill(null(op), xsize, ysize)
    @inbounds for idx in eachindex(xs, ys, zs...)
        x = xs[idx]
        y = ys[idx]
        z = map(z->z[idx], zs)
        xmin ≤ x ≤ xmax || continue
        ymin ≤ y ≤ ymax || continue
        i = clamp(1+floor(Int, xscale*(x-xmin)), 1, xsize)
        j = clamp(1+floor(Int, yscale*(y-ymin)), 1, ysize)
        out[i,j] = update(op, out[i,j], z...)
        nothing
    end
    map(x->value(op, x), out)
end

function _aggregate(c::Canvas, op::AggOp, method::AggThreads, xs, ys, zs...)
    xmin, xmax = xlims(c)
    ymin, ymax = ylims(c)
    xsize, ysize = size(c)
    xmax > xmin || error("require xmax > xmin")
    ymax > ymin || error("require ymax > ymin")
    xscale = xsize / (xmax - xmin)
    yscale = ysize / (ymax - ymin)
    # each thread reduces some of the data separately
    out = fill(null(op), nthreads(), xsize, ysize)
    @threads for idx in eachindex(xs, ys, zs...)
        t = threadid()
        x = @inbounds xs[idx]
        y = @inbounds ys[idx]
        z = map(z->@inbounds(z[idx]), zs)
        xmin ≤ x ≤ xmax || continue
        ymin ≤ y ≤ ymax || continue
        i = clamp(1+floor(Int, xscale*(x-xmin)), 1, xsize)
        j = clamp(1+floor(Int, yscale*(y-ymin)), 1, ysize)
        @inbounds out[t,i,j] = update(op, out[t,i,j], z...)
    end
    # reduce along the thread dimension
    out2 = fill(null(op), xsize, ysize)
    for j in 1:ysize
        for i in 1:xsize
            @inbounds val = out[1,i,j]
            for t in 2:nthreads()
                @inbounds val = merge(op, val, out[t,i,j])
            end
            @inbounds out2[i,j] = val
        end
    end
    map(x->value(op, x), out2)
end

const DEFAULT_SPREAD_MASKS = Matrix{Bool}[
    ones(Int, 1, 1),
    [0 1 0; 1 1 1; 0 1 0],
    [1 1 1; 1 1 1; 1 1 1],
    [0 0 1 0 0; 0 1 1 1 0; 1 1 1 1 1; 0 1 1 1 0; 0 0 1 0 0],
    [0 1 1 1 0; 1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1; 0 1 1 1 0],
    [1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1],
]

function spread(img::AbstractMatrix, r::Integer; masks=DEFAULT_SPREAD_MASKS, opts...)
    Base.require_one_based_indexing(masks)
    spread(img, masks[r+1]; opts...)
end

function spread(img::AbstractMatrix, w::AbstractMatrix; op=max)
    Base.require_one_based_indexing(img, w)
    wsz = size(w)
    all(isodd, wsz) || error("weights must have odd size in each dimension")
    r = map(x->fld(x,2), wsz)
    sz = size(img)
    out = zeros(typeof(zero(eltype(img))*zero(eltype(w))), size(img)...)
    for j in -r[2]:r[2]
        for i in -r[1]:r[1]
            wt = w[r[1]+1+i,r[2]+1+j]
            vout = @view out[i ≥ 0 ? (1+i:end) : (1:end+i), j ≥ 0 ? (j+1:end) : (1:end+j)]
            vimg = @view img[i ≥ 0 ? (1:end-i) : (1-i:end), j ≥ 0 ? (1:end-j) : (1-j:end)]
            vout .= op.(vout, vimg .* Ref(wt))
        end
    end
    out
end

function autospread(img::AbstractMatrix; masks=DEFAULT_SPREAD_MASKS, rmax=length(masks)-1, thresh=0.5, opts...)
    Base.require_one_based_indexing(img, masks)
    0 ≤ rmax < length(masks) || error("rmax out of range")
    n₀ = count(x->!iszero(x), img)
    out = spread(img, masks[1]; opts...)
    for r in 1:rmax
        mask = masks[r+1]
        s = count(x->!iszero(x), mask)
        newout = spread(img, mask; opts...)
        n = count(x->!iszero(x), newout)
        n < thresh * s * n₀ && break
        # out = newout
        # linearly interpolate between out and newout depending on where n is in the
        # interval [thresh * s * n₀, s * n₀]
        p = (n / (s * n₀) - thresh) / (1 - thresh)
        out = @. p * newout + (1 - p) * out
    end
    return out
end

"""
    datashade!(scene, x, y, ...; agg=AggCount(), post=identity, colorize=nothing, opts...)

Shade your data into the given scene at co-ordinates `x` and `y`.

Each pixel is some aggregate of the data points in that pixel, controlled by `agg`. By default, each pixel represents the number of points.

The resulting matrix of aggregated values is then post-processed by `post`, if given. For example you could apply a pointwise map to the values, or use `spread` or `autospread`.

If `colorize` is given, it is a map which converts the post-processed matrix into an image (a matrix of colors). Otherwise, a heatmap is plotted.

All other options are passed on to the plotting command.
"""
function datashade!(ax::Axis, x, y, z...; colorize=nothing, opts...)
    c = lift(ax.finallimits, ax.scene.px_area) do lims, pxarea
        xsize, ysize = Makie.widths(pxarea)
        xmin, ymin = minimum(lims)
        xmax, ymax = maximum(lims)
        Canvas(xmin, xmax, xsize, ymin, ymax, ysize)
    end |> observe_changes
    if isnothing(colorize)
        datashadeheatmap!(ax, c, x, y, z...; opts...)
    else
        datashadeimage!(ax, c, x, y, x...; colorize=colorize, opts...)
    end
end

@recipe(DataShadeHeatmap, canvas, x, y) do scene
    th = default_theme(scene, Heatmap)
    th.agg = AggCount()
    th.post = identity
    th
end

@recipe(DataShadeImage, canvas, x, y) do scene
    th = default_theme(scene, Image)
    th.agg = AggCount()
    th.post = identity
    th.colormap = identity
    th
end

function Makie.plot!(p::DataShadeHeatmap{<:Tuple{Canvas,Vararg}})
    xrange = lift(c -> c.xmin .. c.xmax, p.canvas)
    yrange = lift(c -> c.ymin .. c.ymax, p.canvas)
    pixels = lift(p.agg, p.post, p.converted...) do agg, post, args...
        float(post(aggregate(args...; op=agg)))
    end
    th = Theme()
    for k in keys(default_theme(p.parent, Heatmap))
        th[k] = p[k]
    end
    heatmap!(p, xrange, yrange, pixels; th...)
    p
end

function Makie.plot!(p::DataShadeImage{<:Tuple{Canvas,Vararg}})
    xrange = lift(c -> c.xmin .. c.xmax, p.canvas)
    yrange = lift(c -> c.ymin .. c.ymax, p.canvas)
    pixels = lift(p.agg, p.post, p.colorize, p.converted...) do agg, post, colorize, args...
        colorize(post(aggregate(args...; op=agg)))
    end
    th = Theme()
    for k in keys(default_theme(p.parent, Image))
        th[k] = p[k]
    end
    image!(p, xrange, yrange, pixels; th...)
end

function Makie.data_limits(p::Union{DataShadeHeatmap{<:Tuple{Canvas,Vararg}}, DataShadeImage{<:Tuple{Canvas, Vararg}}})
    xmin, xmax = extrema(p.x[])
    ymin, ymax = extrema(p.y[])
    FRect3D([xmin, ymin, 0], [xmax-xmin, ymax-ymin, 0])
end

end
