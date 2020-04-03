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

data(h::AggOp) = ()
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

struct AggSum{T,D} <: AggOp
    data :: D
end
AggSum(data) = AggSum{eltype(data), typeof(data)}(data)
data(a::AggSum) = (a.data,)
null(a::AggSum{T}) where {T} = zero(T)
embed(a::AggSum{T}, x) where {T} = convert(T, x)
merge(a::AggSum{T}, x::T, y::T) where {T} = x + y
value(a::AggSum{T}, x::T) where {T} = x

struct AggMean{T,N,D} <: AggOp
    data :: D
end
AggMean(data) = AggMean{eltype(data), Int, typeof(data)}(data)
data(a::AggMean) = (a.data,)
null(a::AggMean{T,N}) where {T,N} = zero(T), zero(N)
embed(a::AggMean{T,N}, x) where {T,N} = convert(T,x), oneunit(N)
merge(a::AggMean{T,N}, x::Tuple{T,N}, y::Tuple{T,N}) where {T,N} = x[1]+y[1], x[2]+y[2]
value(a::AggMean{T,N}, x::Tuple{T,N}) where {T,N} = iszero(x[2]) ? oftype(float(x[1])/float(x[2]), NaN) : float(x[1])/float(x[2])

function aggregate(c::Canvas, xs, ys; op::AggOp=AggCount())
    xmin, xmax = xlims(c)
    ymin, ymax = ylims(c)
    xsize, ysize = size(c)
    xmax > xmin || error("require xmax > xmin")
    ymax > ymin || error("require ymax > ymin")
    xscale = 1 / (xmax - xmin)
    yscale = 1 / (ymax - ymin)
    out = fill(null(op), xsize, ysize)
    for row in zip(xs, ys, data(op)...)
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

canvas_node(c::Node{<:Canvas}) = c
canvas_node(c::Canvas) = Node{Canvas}(c)

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

set_limits!(scene, how, xs, ys) =
    if (how === :extrema) || (how === true)
        xlims!(scene, extrema(xs))
        ylims!(scene, extrema(ys))
        return
    elseif (how === nothing) || (how === false)
        return
    elseif how isa Node
        scene.limits = how
    else
        scene.limits[] = how
    end

function datashade!(scene::AbstractScene, xs, ys; op::Union{<:AggOp,Node{<:AggOp}}=AggCount(), post=identity, limits=true, xautolimits=false, yautolimits=false, opts...)
    op = op isa Node ? op : Node{AggOp}(op)
    post = post isa Node ? post : Node{Function}(post)
    c = canvas_node(scene)
    xrange = lift(c -> c.xmin .. c.xmax, c)
    yrange = lift(c -> c.ymin .. c.ymax, c)
    pixels = lift(c, op, post) do c, op, post
        float(post(aggregate(c, xs, ys; op=op)))
    end
    plot = heatmap!(scene, xrange, yrange, pixels; xautolimits=xautolimits, yautolimits=yautolimits, opts...)
    set_limits!(scene, limits, xs, ys)
    nodes = (op=op, canvas=c, xrange=xrange, yrange=yrange, pixels=pixels, post=post)
    plot, nodes
end

end # module
