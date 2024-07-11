using Makie
using GLMakie
using CairoMakie
CairoMakie.activate!(type = "svg")
using SwarmMakie
using DataFrames, CSV
using LsqFit
using StatsBase

function figS3A()
    data = DataFrame(CSV.File("/mnt/b/Papers/FAP dispersal/Data/Figure S3/Figure S3A/Compiled.csv"))
	@. model(x, p) = p[1] + (x^p[4])*(p[2]-p[1])/(x^p[4]+p[3]^p[4]) 
    fig = Figure(size=(10*72, 3.5*72), fonts = (; bold = "TeX Gyre Heros Makie"))

    x = data[!, "MgE Conc"]
    y = data[!, "MgE Avg"]
    p0 = [0.0, maximum(y), 0.01, 2]
    lb = [0.0, 0.0, 0.0, 0.0]
    ub = [Inf, Inf, maximum(x), 20.0]
    fit = curve_fit(model, x, y, p0, lower=lb, upper=ub)
    pstar = coef(fit)
    sigma = stderror(fit)
	xbase = collect(range(minimum(x), maximum(x), 100))
    ax = Axis(fig[1,1])
    scatter!(ax, data[!, "MgE Conc"], data[!, "MgE Avg"], strokecolor=:black, strokewidth=1, color=:white)
    errorbars!(ax, data[!, "MgE Conc"], data[!, "MgE Avg"], data[!, "MgE SD"], color=:black)
	lines!(ax, xbase, model.(xbase, (pstar,)), color=:black)
    text!(ax, 0.1, 0, text = "Ec50 = \n"*string(round(pstar[3], digits=2))*"±"*string(round(sigma[3], digits=2))*" μM", 
          align = (:left, :bottom), color=:black)
    ax.rightspinevisible=false
    ax.topspinevisible=false
    ax.xgridvisible=false
    ax.ygridvisible=false
    ax.ylabel = "Fluorescence (a.u.)"
    ax.xlabel = "Concentration (μM)"
    ax.title = "MGe"
    
    x = data[!, "Mg2P Conc"]
    y = data[!, "Mg2P Avg"]
    p0 = [0.0, maximum(y), 0.01, 2]
    lb = [0.0, 0.0, 0.0, 0.0]
    ub = [Inf, Inf, maximum(x), 20.0]
    fit = curve_fit(model, x, y, p0, lower=lb, upper=ub)
    pstar = coef(fit)
    sigma = stderror(fit)
	xbase = collect(range(minimum(x), maximum(x), 100))
    ax2 = Axis(fig[1,2])
    scatter!(ax2, data[!, "Mg2P Conc"], data[!, "Mg2P Avg"], strokecolor=:black, strokewidth=1, color=:white)
    errorbars!(ax2, data[!, "Mg2P Conc"], data[!, "Mg2P Avg"], data[!, "Mg2P SD"], color=:black)
	lines!(ax2, xbase, model.(xbase, (pstar,)), color=:black)
    text!(ax2, 0.1, 0, text = "Ec50 = \n"*string(round(pstar[3], digits=2))*"±"*string(round(sigma[3], digits=2))*" μM", 
          align = (:left, :bottom), color=:black)
    ax2.rightspinevisible=false
    ax2.topspinevisible=false
    ax2.xgridvisible=false
    ax2.ygridvisible=false
    ax2.ylabel = ""
    ax2.xlabel = "Concentration (μM)"
    ax2.title = "MG-2P"
    
    x = data[!, "MHN Conc"]
    y = data[!, "MHN Avg"]
    p0 = [0.0, maximum(y), 0.01, 2]
    lb = [0.0, 0.0, 0.0, 0.0]
    ub = [Inf, Inf, maximum(x), 20.0]
    fit = curve_fit(model, x, y, p0, lower=lb, upper=ub)
    pstar = coef(fit)
    sigma = stderror(fit)
	xbase = collect(range(minimum(x), maximum(x), 100))
    ax2 = Axis(fig[1,3])
    scatter!(ax2, data[!, "MHN Conc"], data[!, "MHN Avg"], strokecolor=:black, strokewidth=1, color=:white)
    errorbars!(ax2, data[!, "MHN Conc"], data[!, "MHN Avg"], data[!, "MHN SD"], color=:black)
	lines!(ax2, xbase, model.(xbase, (pstar,)), color=:black)
    text!(ax2, 1, 0, text = "Ec50 = \n"*string(round(pstar[3], digits=2))*"±"*string(round(sigma[3], digits=2))*" μM", 
          align = (:left, :bottom), color=:black)
    ax2.rightspinevisible=false
    ax2.topspinevisible=false
    ax2.xgridvisible=false
    ax2.ygridvisible=false
    ax2.ylabel = ""
    ax2.xlabel = "Concentration (μM)"
    ax2.title = "MHNe"
    save("/mnt/b/Papers/FAP dispersal/Figures/FigS3A.svg", fig)
end

function figS3B()
    data = DataFrame(CSV.File("/mnt/b/Papers/FAP dispersal/Data/Figure S3/Figure S3B/FigS3B.csv"))
    averages = Float64.([mean(data[!, i]) for i in 1:3]) ./ mean(data[!, 1])
    mins = Float64.([minimum(data[!, i]) for i in 1:3]) ./ mean(data[!, 1])
    maxes = Float64.([maximum(data[!, i]) for i in 1:3]) ./ mean(data[!, 1])
    conditions = ["Background", "MGe", "MG-2P"]
    fig = Figure(size=(6*72, 3.5*72))
    ax = Axis(fig[1,1])
    ax2 = Axis(fig[1,2])
	category_num_swarm = Int.(repeat(1:3, inner = 3))
	category_num = Int.(1:3)
	category_num_swarm = Int.(repeat(1:3, inner=3))
	colormap1 = [[:black]; [:magenta, :magenta]]
	colormap2 = [[:white]; [:magenta, :magenta]]
    data1 = vcat([Float64.(data[!, i]) for i in 1:3]...) ./ mean(data[!, 1])
	crossbar!(ax, category_num, averages, mins, maxes; color=:white, midlinecolor=colormap1, colormap1, colorrange=(1,3))
	plt = beeswarm!(ax, category_num_swarm, data1, color = category_num_swarm, algorithm=UniformJitter(), strokewidth=1)
	plt.colormap[] = colormap2
	ax.xticks=(1:3, conditions)
	ax.xticklabelrotation=45
	ax.xlabel=""
	ax.ylabel="Relative fluorescence (a.u.)"
	ax.title=""
	ax.rightspinevisible = false
	ax.topspinevisible = false
	ax.xgridvisible = false
	ax.ygridvisible = false
	
    averages = Float64.([mean(data[!, i]) for i in 4:5]) ./ mean(data[!, 4]) 
    mins = Float64.([minimum(data[!, i]) for i in 4:5]) ./ mean(data[!, 4])
    maxes = Float64.([maximum(data[!, i]) for i in 4:5]) ./ mean(data[!, 4])
    conditions = ["Background", "MHNe"]
	category_num_swarm = Int.(repeat(1:2, inner = 3))
	category_num = Int.(1:2)
	category_num_swarm = Int.(repeat(1:2, inner=3))
	colormap1 = [[:black]; [:cyan2]]
	colormap2 = [[:white]; [:cyan2]]
    data2 = vcat([Float64.(data[!, i]) for i in 4:5]...) ./ mean(data[!, 4])
	crossbar!(ax2, category_num, averages, mins, maxes; color=:white, midlinecolor=colormap1, colormap1, colorrange=(1,2))
	plt = beeswarm!(ax2, category_num_swarm, data2, color = category_num_swarm, algorithm=UniformJitter(), strokewidth=1)
	plt.colormap[] = colormap2
	ax2.xticks=(1:2, conditions)
	ax2.xticklabelrotation=45
	ax2.xlabel=""
	ax2.ylabel=""
	ax2.title=""
	ax2.rightspinevisible = false
	ax2.topspinevisible = false
	ax2.xgridvisible = false
	ax2.ygridvisible = false
    ylims!(ax2, 0, nothing)
    save("/mnt/b/Papers/FAP dispersal/Figures/FigS3B.svg", fig)
end

function figS3C()
    data_mg2p = DataFrame(CSV.File("/mnt/b/Papers/FAP dispersal/Data/Figure S3/Figure S3C/Add fluorogens/Results/Normalized.csv"))
    data_mge = DataFrame(CSV.File("/mnt/b/Papers/FAP dispersal/Data/Figure S3/Figure S3C/Add fluorogens/Results/Normalized_MGe.csv"))
    fig = Figure(size=(5*72, 3.5*72))
    ax = Axis(fig[1,1])
    lines!(ax, data_mg2p[!, "Time"], data_mg2p[!, "AVG"], color=:magenta)
    band!(ax, data_mg2p[!, "Time"], data_mg2p[!, "AVG"] - data_mg2p[!, "STDEV"], data_mg2p[!, "AVG"] + data_mg2p[!, "STDEV"],
          color=(:magenta,0.3))
    scatter!(ax, data_mg2p[!, "Time"], data_mg2p[!, "AVG"], color=:white, strokecolor=:black, strokewidth=1, label="MG-2P")
    lines!(ax, data_mge[!, "Time"], data_mge[!, "AVG"], color=:magenta)
    band!(ax, data_mge[!, "Time"], data_mge[!, "AVG"] - data_mge[!, "STDEV"], data_mge[!, "AVG"] + data_mge[!, "STDEV"],
          color=(:magenta,0.3))
    scatter!(ax, data_mge[!, "Time"], data_mge[!, "AVG"], color=:magenta, strokecolor=:black, strokewidth=1, label="MGe")
    ax.rightspinevisible=false
    ax.topspinevisible=false
    ax.xgridvisible=false
    ax.ygridvisible=false
    ax.ylabel = "Normalized fluorescence (a.u.)"
    ax.xlabel = "Time (min)"
    fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)
    save("/mnt/b/Papers/FAP dispersal/Figures/FigS3C.svg", fig)
end

figS3C()
