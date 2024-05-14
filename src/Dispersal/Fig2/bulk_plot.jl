using Makie
using GLMakie
using CairoMakie
CairoMakie.activate!(type="svg")
using LaTeXStrings
using StatsBase
using DelimitedFiles
using Makie.Colors

function main()
    plots_folder = "/mnt/h/Dispersal/Plots"
    files = [f for f in readdir(plots_folder, join=true) if occursin("processed.csv", f)]
    plot_xlabel = "Time (h)"
    plot_ylabel = "Biovolume (a.u.)" 
    fig=Figure(size=(5*72, 2.5*72))
    ax = Axis(fig[1,1])
    WT_seen = false
    cheY_seen = false
    lapG_seen = false
    rbmB_seen = false
    for file in files
        data = readdlm(file, ',', Int)[1:end,1]
        first_idx = argmax(data)
        end_idx = min(first_idx + 45, length(data))
        data_subset = data[first_idx:end_idx] 
        data = data_subset
        data_norm =  data ./ data[1]
        data = data_norm
        lab = split(basename(file), "_")[1]
        if lab == "WT"
            color = 1 
            colormap = :Pastel1_4
            if WT_seen
                condition = nothing 
            else
                WT_seen = true
                condition = "Wild-type"
            end
        elseif lab == "cheY"
            color = 2
            colormap = :Pastel1_4
            if cheY_seen
                condition = nothing 
            else
                cheY_seen = true
                condition = rich("Δ", rich("cheY"; font=:italic))
            end
        elseif lab == "rbmB"
            color = 3
            colormap = :Pastel1_4
            if rbmB_seen
                condition = nothing 
            else
                rbmB_seen = true
                condition = rich("Δ", rich("rbmB"; font=:italic))
            end
        elseif lab == "lapG"
            color = 4 
            colormap = :Pastel1_4
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
        lines!(ax, xaxis, data, color=color, colormap=colormap, colorrange = (1, 4), label=condition)
    end
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
end

main()
