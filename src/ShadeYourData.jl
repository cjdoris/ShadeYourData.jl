module ShadeYourData

export datashade!

using AbstractPlotting, StatsBase, Requires

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

abstract type AggOp end

update(a::AggOp, x, args...) = merge(a, x, embed(a, args...))

struct AggCount{T} <: AggOp end
AggCount() = AggCount{Int}()
null(a::AggCount{T}) where {T} = zero(T)
embed(a::AggCount{T}) where {T} = oneunit(T)
merge(a::AggCount{T}, x::T, y::T) where {T} = x + y
value(a::AggCount{T}, x::T) where {T} = x

struct AggAny{T} <: AggOp end
AggAny() = AggAny{Bool}()
null(a::AggAny{T}) where {T} = zero(T)
embed(a::AggAny{T}) where {T} = oneunit(T)
merge(a::AggAny{T}, x::T, y::T) where {T} = max(x, y)
value(a::AggAny{T}, x::T) where {T} = x

struct AggSum{T} <: AggOp end
AggSum() = AggSum{Float64}()
null(a::AggSum{T}) where {T} = zero(T)
embed(a::AggSum{T}, x) where {T} = convert(T, x)
merge(a::AggSum{T}, x::T, y::T) where {T} = x + y
value(a::AggSum{T}, x::T) where {T} = x

struct AggMean{T,N} <: AggOp end
AggMean{T}() where {T} = AggMean{T,Int}()
AggMean() = AggMean{Float64}()
null(a::AggMean{T,N}) where {T,N} = zero(T), zero(N)
embed(a::AggMean{T,N}, x) where {T,N} = convert(T,x), oneunit(N)
merge(a::AggMean{T,N}, x::Tuple{T,N}, y::Tuple{T,N}) where {T,N} = x[1]+y[1], x[2]+y[2]
value(a::AggMean{T,N}, x::Tuple{T,N}) where {T,N} = float(x[1]) / float(x[2])

function aggregate(c::Canvas, xs, ys, zs...; op::AggOp=AggCount())
    xmin, xmax = xlims(c)
    ymin, ymax = ylims(c)
    xsize, ysize = size(c)
    xmax > xmin || error("require xmax > xmin")
    ymax > ymin || error("require ymax > ymin")
    xscale = 1 / (xmax - xmin)
    yscale = 1 / (ymax - ymin)
    out = fill(null(op), xsize, ysize)
    for row in zip(xs, ys, zs...)
        x = row[1]
        y = row[2]
        z = row[3:end]
        xmin ≤ x ≤ xmax || continue
        ymin ≤ y ≤ ymax || continue
        i = clamp(1+floor(Int, xsize*xscale*(x-xmin)), 1, xsize)
        j = clamp(1+floor(Int, ysize*yscale*(y-ymin)), 1, ysize)
        @inbounds out[i,j] = update(op, out[i,j], z...)
    end
    map(x->value(op, x), out)
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
        n ≥ thresh * s * n₀ || return out
        out = newout
    end
    return out
end

function autospread_old(img::AbstractMatrix, rmax::Integer=5, op=max; thresh=0.2)
    outmax = spread(img, rmax, op)
    nummax = count(x->!iszero(x), outmax)
    for r in 0:rmax-1
        out = spread(img, r, op)
        num = count(x->!iszero(x), out)
        if num ≥ thresh * nummax
            return out
        end
    end
    return outmax
end

canvas_node(c::Canvas) = Node{Canvas}(c)
canvas_node(c::Node{<:Canvas}) = c

function __init__()
    @require MakieLayout="5a521ce4-ebb9-4793-b5b7-b334dfe8393c" begin
        canvas_node(ax::MakieLayout.LAxis) = lift(ax.limits, ax.scene.px_area) do lims, pxarea
            xsize, ysize = widths(pxarea)
            xmin, ymin = minimum(lims)
            xmax, ymax = maximum(lims)
            Canvas(xmin, xmax, xsize, ymin, ymax, ysize)
        end
    end    
end

"""
    datashade!(scene, x, y, ...; agg=AggCount(), post=identity, colorize=nothing, opts...)

Shade your data into the given scene at co-ordinates `x` and `y`.

Each pixel is some aggregate of the data points in that pixel, controlled by `agg`. By default, each pixel represents the number of points.

The resulting matrix of aggregated values is then post-processed by `post`, if given. For example you could apply a pointwise map to the values, or use `spread` or `autospread`.

If `colorize` is given, it is a map which converts the post-processed matrix into an image (a matrix of colors). Otherwise, a heatmap is plotted.

All other options are passed on to the plotting command.
"""
function datashade!(scene::AbstractScene, x, y, z...; colorize=nothing, opts...)
    c = canvas_node(scene)
    if isnothing(colorize)
        datashadeheatmap!(scene, c, x, y, z...; opts...)
    else
        datashadeimage!(scene, c, x, y, x...; colorize=colorize, opts...)
    end
end

@recipe(DataShadeHeatmap, canvas, x, y, z) do scene
    th = default_theme(scene, Heatmap)
    th.agg = AggCount()
    th.post = identity
    th
end

@recipe(DataShadeImage, canvas, x, y, z) do scene
    th = default_theme(scene, Image)
    th.agg = AggCount()
    th.post = identity
    th.colormap = identity
    th
end

function AbstractPlotting.plot!(p::DataShadeHeatmap{<:Tuple{Canvas,Vararg}})
    global LAST_PLOT = p
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
end

function AbstractPlotting.plot!(p::DataShadeImage{<:Tuple{Canvas,Vararg}})
    global LAST_PLOT = p
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

function AbstractPlotting.data_limits(p::Union{DataShadeHeatmap{<:Tuple{Canvas,Vararg}}, DataShadeImage{<:Tuple{Canvas, Vararg}}})
    xmin, xmax = extrema(p.x[])
    ymin, ymax = extrema(p.y[])
    FRect3D([xmin, ymin, 0], [xmax-xmin, ymax-ymin, 0])
end

end # module
