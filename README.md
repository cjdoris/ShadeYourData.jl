# ShadeYourData.jl

**ShadeYourData.jl** *n.*
1. A cheap clone of [datashader.org](https://datashader.org) for Julia.
2. Interactive plotting of millions of data points.

## Getting started

Install like so: `] add https://github.com/cjdoris/ShadeYourData.jl`.

See the example below, or read the docstring for `datashade!`.

## Example

In this example, we plot all 2.6M UK postcodes from [this dataset](https://www.doogal.co.uk/PostcodeDownloads.php).

![ShadeYourData.jl example](https://raw.githubusercontent.com/cjdoris/ShadeYourData.jl/master/example.png)

You can watch it in action [here](http://www.youtube.com/watch?v=svG6fCjVbEg "ShadeYourData.jl demo"). Notice that when zoomed out, the density of points is shown via a colormap, but as we zoom in and the points get separated, we plot them bigger so you can see them individually.

```julia
# ShadeYourData works best with MakieLayout
using CSV, DataFrames, Makie, MakieLayout, ShadeYourData

# load the dataset
df = CSV.read("/path/to/postcodes.csv", copycols=true)
dropmissing!(df, [:Easting, :Northing])

# set up a new scene
scene, layout = layoutscene()
ax = layout[1,1] = LAxis(scene)

# `datashade!` is the main plotting routine we provide
ds = datashade!(ax, df.Easting, df.Northing, colormap=:dense)

# post-process the aggregated image by applying an autospread (makes points bigger when they are spread out) and rescaling by `log1p` to lessen the extremal values
ds.post = x -> log1p.(ShadeYourData.autospread(x))

# force the data aspect ratio to be 1
# you can achieve the same with just `ax.autolimitaspect = true`, but currently zooming doesn't work well when you do this
ax.aspect = DataAspect()
xlims!(ax, -300_000, 900_000)
ylims!(ax, -100_000, 1_100_000)

# take a look
display(scene)
```
