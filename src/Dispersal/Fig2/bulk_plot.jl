using Makie
using GLMakie
using CairoMakie
CairoMakie.activate!(type="svg")
using LaTeXStrings
using StatsBase
using DelimitedFiles
using Makie.Colors
using SwarmMakie

function main()
    plots_folder = "/mnt/h/Dispersal/Plots"
    files = [f for f in readdir(plots_folder, join=true) if occursin("processed.csv", f) && !occursin("rbmA", f)]
    plot_xlabel = "Time (h)"
    plot_ylabel = "Biovolume (a.u.)" 
    fig=Figure(size=(5*72, 2.5*72))
    ax = Axis(fig[1,1])
    WT_seen = false
    cheY_seen = false
    lapG_seen = false
    rbmB_seen = false
    rbmA_files = [f for f in readdir(plots_folder, join=true) if occursin("rbmA", f) && occursin("processed.csv", f)]
    WT = []
    cheY = []
    lapG = []
    rbmB = []
    rbmA = []
    WTs = []
    cheYs = []
    lapGs = []
    rbmBs = []
    rbmAs = []
	column_names = ["WT", "cheY", "lapG", "rbmB", "rbmA"]
    for file in rbmA_files
        data = readdlm(file, ',', Int)[1:end,1]
        peak_biofilm = maximum(data)
        push!(rbmA, peak_biofilm)
		data = readdlm(file, ',', Int)[1:end,1] .+ 0.01
		data_max = maximum(data)
		fc = data ./ data_max
		min_fc = minimum(fc)
		fc = fc .- min_fc
		end_idx = findfirst(x->x<0.001, fc)
		first_idx = end_idx - 45
        data_subset = data[first_idx:end_idx] 
        data = data_subset
        data_norm =  data ./ data[1]
        data = data_norm
    end
    for file in files
        data = readdlm(file, ',', Int)[1:end,1]
        first_idx = argmax(data)
        peak_biofilm = maximum(data)
        end_idx = min(first_idx + 45, length(data))
        data_subset = data[first_idx:end_idx] 
        data = data_subset
        data_norm =  data ./ data[1]
        data = data_norm
        lab = split(basename(file), "_")[1]
        if lab == "WT"
            push!(WT, peak_biofilm)
            color = 1 
            colormap = [[:black];  Makie.wong_colors()[1:3]]
            if WT_seen
                condition = nothing 
            else
                WT_seen = true
                condition = "Wild-type"
            end
        elseif lab == "cheY"
            push!(cheY, peak_biofilm)
            color = 2
            colormap = [[:black];  Makie.wong_colors()[1:3]]
            if cheY_seen
                condition = nothing 
            else
                cheY_seen = true
                condition = rich("Δ", rich("cheY"; font=:italic))
            end
        elseif lab == "rbmB"
            push!(rbmB, peak_biofilm)
            color = 4
            colormap = [[:black];  Makie.wong_colors()[1:3]]
            if rbmB_seen
                condition = nothing 
            else
                rbmB_seen = true
                condition = rich("Δ", rich("rbmB"; font=:italic))
            end
        elseif lab == "lapG"
            push!(lapG, peak_biofilm)
            color = 3 
            colormap = [[:black];  Makie.wong_colors()[1:3]]
            if lapG_seen
                condition = nothing 
            else
                lapG_seen = true
                condition = rich("Δ", rich("lapG"; font=:italic))
            end
        end
        xs = 0:6:length(data)-1
        xaxis = 1:length(data)
        xaxis /= 6
        lines!(ax, xaxis, data, color=color, colormap=colormap, colorrange = (1,4), label=condition)
    end
    conditions = ["Wild-type", rich("Δ", rich("cheY"; font=:italic)), rich("Δ", rich("lapG"; font=:italic)), 
                  rich("Δ", rich("rbmB"; font=:italic)), rich("Δ", rich("rbmA"; font=:italic))]
    fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)
    ax.xlabel = plot_xlabel
    ax.ylabel = plot_ylabel
    ax.title = ""
    ax.rightspinevisible = false
    ax.topspinevisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    ylims!(ax, 0, nothing)
    save(plots_folder*"/WT_cheY_lapG_rbmB_bulk.svg", fig)

	data = vcat(WT, cheY, lapG, rbmB, rbmA) .* 0.065^3
    averages = [mean(data[1:5]), mean(data[6:10]), mean(data[11:15]), mean(data[16:20]), mean(data[21:25])]
    maxes = [maximum(data[1:5]), maximum(data[6:10]), maximum(data[11:15]), maximum(data[16:20]), maximum(data[21:25])]
    mins = [minimum(data[1:5]), minimum(data[6:10]), minimum(data[11:15]), minimum(data[16:20]), minimum(data[21:25])]
    category_num_swarm = Int.(repeat(1:5, inner = 5))
    fig = Figure(size=(4*72, 3*72))
	category_num = Int.(1:5)
	category_num_swarm = Int.(repeat(1:5, inner=5))
	ax = Axis(fig[1, 1])
	crossbar!(ax, category_num, averages, mins, maxes; 
			  color=:white, midlinecolor=:black)
    plt = beeswarm!(ax, category_num_swarm, Float64.(data), color = :white, algorithm=UniformJitter(), strokewidth=1)
    ax.xticks=(1:5, conditions)
    ax.xticklabelrotation=45
    ax.xlabel=""
    ax.ylabel=rich("Peak biofilm volume \n (μm", superscript("3"), ")")
    ax.title=""
    ax.rightspinevisible = false
    ax.topspinevisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    ylims!(ax, 0, nothing)
    save(plots_folder*"/peak_biofilm_size.svg", fig)
end

main()
