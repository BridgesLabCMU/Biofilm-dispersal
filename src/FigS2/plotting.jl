using Makie
using GLMakie
using CairoMakie
CairoMakie.activate!(type = "svg")
using DataFrames, CSV
using StatsBase

function figS2A()
    dataS2A = "/mnt/b/Papers/FAP dispersal/Data/Figure S2/Figure S2A/GrowthData.csv"
    dataS2A = DataFrame(CSV.File(dataS2A)) 
    fig = Figure(size = (5*72,3*72))
    ax = Axis(fig[1,1])
    lines!(ax, dataS2A[!, "Time"], dataS2A[!, "AVG_310"], color=:black)
    lines!(ax, dataS2A[!, "Time"], dataS2A[!, "AVG_1548"], color=:black)
    lines!(ax, dataS2A[!, "Time"], dataS2A[!, "AVG_1717"], color=:black)
    band!(ax, dataS2A[!, "Time"], dataS2A[!, "AVG_310"]-dataS2A[!, "STDEV_310"], dataS2A[!, "AVG_310"]+dataS2A[!, "STDEV_310"], color=(:black, 0.2))
    band!(ax, dataS2A[!, "Time"], dataS2A[!, "AVG_1548"]-dataS2A[!, "STDEV_1548"], dataS2A[!, "AVG_1548"]+dataS2A[!, "STDEV_1548"], color=(:black, 0.2))
    band!(ax, dataS2A[!, "Time"], dataS2A[!, "AVG_1717"]-dataS2A[!, "STDEV_1717"], dataS2A[!, "AVG_1717"]+dataS2A[!, "STDEV_1717"], color=(:black, 0.2))
    scatter!(ax, dataS2A[!, "Time"], dataS2A[!, "AVG_310"], color=:white, marker=:circle,  strokewidth=1, label=rich("mNG"; font=:italic))
    scatter!(ax, dataS2A[!, "Time"], dataS2A[!, "AVG_1548"], color=:white, marker=:rect,  strokewidth=1, label=rich("mNG(Y69G)-dL5"; font=:italic))
    scatter!(ax, dataS2A[!, "Time"], dataS2A[!, "AVG_1717"], color=:white, marker=:utriangle, strokewidth=1, label=rich("SS-dL5"; font=:italic))
    fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)
    ax.xlabel = "Time (h)"
    ax.ylabel = rich("OD", subscript("600"))
	ax.rightspinevisible = false
	ax.topspinevisible = false
	ax.xgridvisible = false
	ax.ygridvisible = false
    save("/mnt/b/Papers/FAP dispersal/Figures/FigS2A.svg", fig)
end

function figS2B()
    dataS2B = "/mnt/b/Papers/FAP dispersal/Data/Figure S2/Figure S2B/GrowthData.csv"
    dataS2B = DataFrame(CSV.File(dataS2B))
    fig = Figure(size = (5*72,3*72))
    ax = Axis(fig[1,1])
    lines!(ax, dataS2B[!, "Time"], dataS2B[!, "AVG_1717_MgE_C"], color=:black)
    lines!(ax, dataS2B[!, "Time"], dataS2B[!, "AVG_1717_MgE"], color=:magenta)
    band!(ax, dataS2B[!, "Time"], dataS2B[!, "AVG_1717_MgE_C"]-dataS2B[!, "STDEV_1717_MgE_C"], dataS2B[!, "AVG_1717_MgE_C"]+dataS2B[!, "STDEV_1717_MgE_C"], color=(:black, 0.2))
    band!(ax, dataS2B[!, "Time"], dataS2B[!, "AVG_1717_MgE"]-dataS2B[!, "STDEV_1717_MgE"], dataS2B[!, "AVG_1717_MgE"]+dataS2B[!, "STDEV_1717_MgE"], color=(:magenta, 0.2))
    scatter!(ax, dataS2B[!, "Time"], dataS2B[!, "AVG_1717_MgE_C"], color=:black, strokewidth=1, label="No fluorogen")
    scatter!(ax, dataS2B[!, "Time"], dataS2B[!, "AVG_1717_MgE"], color=:magenta, strokewidth=1, label="1 µM MGe")
    fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)
    ax.xlabel = "Time (h)"
    ax.ylabel = rich("OD", subscript("600"))
	ax.rightspinevisible = false
	ax.topspinevisible = false
	ax.xgridvisible = false
	ax.ygridvisible = false
    save("/mnt/b/Papers/FAP dispersal/Figures/FigS2B1.svg", fig)
    
    fig = Figure(size = (5*72,3*72))
    ax = Axis(fig[1,1])
    lines!(ax, dataS2B[!, "Time"], dataS2B[!, "AVG_1717_MHN_C"], color=:black)
    lines!(ax, dataS2B[!, "Time"], dataS2B[!, "AVG_1717_MHN"], color=:cyan2)
    band!(ax, dataS2B[!, "Time"], dataS2B[!, "AVG_1717_MHN_C"]-dataS2B[!, "STDEV_1717_MHN_C"], dataS2B[!, "AVG_1717_MHN_C"]+dataS2B[!, "STDEV_1717_MHN_C"], color=(:black, 0.2))
    band!(ax, dataS2B[!, "Time"], dataS2B[!, "AVG_1717_MHN"]-dataS2B[!, "STDEV_1717_MHN"], dataS2B[!, "AVG_1717_MHN"]+dataS2B[!, "STDEV_1717_MHN"], color=(:cyan2, 0.2))
    scatter!(ax, dataS2B[!, "Time"], dataS2B[!, "AVG_1717_MHN_C"], color=:black, strokewidth=1, label="No fluorogen")
    scatter!(ax, dataS2B[!, "Time"], dataS2B[!, "AVG_1717_MHN"], color=:cyan2, strokewidth=1, label="5 µM MHNe")
    fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)
    ax.xlabel = "Time (h)"
    ax.ylabel = rich("OD", subscript("600"))
	ax.rightspinevisible = false
	ax.topspinevisible = false
	ax.xgridvisible = false
	ax.ygridvisible = false
    save("/mnt/b/Papers/FAP dispersal/Figures/FigS2B2.svg", fig)
    
    fig = Figure(size = (5*72,3*72))
    ax = Axis(fig[1,1])
    lines!(ax, dataS2B[!, "Time"], dataS2B[!, "AVG_1548"], color=:black)
    lines!(ax, dataS2B[!, "Time"], dataS2B[!, "AVG_1548_MG2P"], color=:magenta)
    band!(ax, dataS2B[!, "Time"], dataS2B[!, "AVG_1548"]-dataS2B[!, "STDEV_1548"], dataS2B[!, "AVG_1548"]+dataS2B[!, "STDEV_1548"], color=(:black, 0.2))
    band!(ax, dataS2B[!, "Time"], dataS2B[!, "AVG_1548_MG2P"]-dataS2B[!, "STDEV_1548_MG2P"], dataS2B[!, "AVG_1548_MG2P"]+dataS2B[!, "STDEV_1548_MG2P"], color=(:magenta, 0.2))
    scatter!(ax, dataS2B[!, "Time"], dataS2B[!, "AVG_1548"], color=:black, strokewidth=1, label="No fluorogen")
    scatter!(ax, dataS2B[!, "Time"], dataS2B[!, "AVG_1548_MG2P"], color=:white, strokewidth=1, label="1 µM MG-2P")
    fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)
    ax.xlabel = "Time (h)"
    ax.ylabel = rich("OD", subscript("600"))
	ax.rightspinevisible = false
	ax.topspinevisible = false
	ax.xgridvisible = false
	ax.ygridvisible = false
    save("/mnt/b/Papers/FAP dispersal/Figures/FigS2B3.svg", fig)
end

figS2B()
