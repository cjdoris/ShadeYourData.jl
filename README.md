# ShadeYourData.jl

**ShadeYourData.jl** *n.*
1. A cheap clone of [datashader.org](https://datashader.org) for Julia.
2. Interactive plotting of millions of data points.

## Getting started

Install like so: `] add https://github.com/cjdoris/ShadeYourData.jl`.

See the example below, or read the docstring for `datashade!`.

## Example

In this example, we plot all 2.6M UK postcodes from [this dataset](https://www.doogal.co.uk/PostcodeDownloads.php).

Click the following link to see it in action. Notice that when zoomed out, the density of points is shown via a colormap, but as we zoom in and the points get separated, we plot them bigger so you can see them individually.

[![ShadeYourData.jl demo thumbnail](http://img.youtube.com/vi/svG6fCjVbEg/0.jpg)](http://www.youtube.com/watch?v=svG6fCjVbEg "ShadeYourData.jl demo")

```julia
# ShadeYourData works best with MakieLayout
using CSV, Makie, MakieLayout, ShadeYourData

# load the dataset
df = CSV.read("/path/to/postcodes.csv", copycols=true)
dropmissing!(df, [:Easting, :Northing])

# set up a new scene
scene, layout = layoutscene()
ax = layout[1,1] = LAxis(scene, aspect=DataAspect())

# `datashade!` is the main plotting routine we provide, it returns the plot and a namedtuple of `Node`s that the plot depends on
ds, dsattr = datashade!(ax, df.Easting, df.Northing, colormap=:dense, limits=false)

# manually set the limits
ax.targetlimits[] = FRect(-300_000, -100_000, 1_200_000, 1_200_000)

# post-process the aggregated image by applying an autospread (makes points bigger when they are spread out) and rescaling by `log1p` to lessen the extremal values
dsattr.post[] = x -> log1p.(ShadeYourData.autospread(x))

# take a look
display(scene)
```
