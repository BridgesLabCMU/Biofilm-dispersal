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
    files = [f for f in readdir(plots_folder, join=true) if occursin("processed.csv", f) && occursin("WT", f)]
    plot_xlabel = "Time (h)"
    plot_ylabel = "Biovolume (a.u.)" 
    fig=Figure(size=(3.5*72, 2.5*72))
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
            colormap = [[:black];  Makie.wong_colors()[1:3]]
            if WT_seen
                condition = nothing 
            else
                WT_seen = true
                condition = "Wild-type"
            end
        elseif lab == "cheY"
            color = 2
            colormap = [[:black];  Makie.wong_colors()[1:3]]
            if cheY_seen
                condition = nothing 
            else
                cheY_seen = true
                condition = rich("Δ", rich("cheY"; font=:italic))
            end
        elseif lab == "rbmB"
            color = 4 
            colormap = [[:black];  Makie.wong_colors()[1:3]]
            if rbmB_seen
                condition = nothing 
            else
                rbmB_seen = true
                condition = rich("Δ", rich("rbmB"; font=:italic))
            end
        elseif lab == "lapG"
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
        if occursin("replicate2", file) || occursin("replicate1", file)
            lines!(ax, xaxis, data, color=:black, label=condition)
        else
            lines!(ax, xaxis, data, color=color, colormap=colormap, colorrange = (1, 6), label=condition)
        end
    end
    ax.xlabel = plot_xlabel
    ax.ylabel = plot_ylabel
    ax.title = ""
    ax.rightspinevisible = false
    ax.topspinevisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    ylims!(ax, 0, nothing)
    save(plots_folder*"/WT_bulk_colored.svg", fig)
end

main()
