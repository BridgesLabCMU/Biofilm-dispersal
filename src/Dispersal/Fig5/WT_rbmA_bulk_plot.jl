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
    files = [f for f in readdir(plots_folder, join=true) if occursin("processed.csv", f) && (occursin("WT", f) || occursin("rbmA", f))]
    plot_xlabel = "Time (h)"
    plot_ylabel = "Biovolume (a.u.)" 
    fig=Figure(size=(5*72, 2.5*72))
    ax = Axis(fig[1,1])
    WT_seen = false
    rbmA_seen = false
    for file in files
        lab = split(basename(file), "_")[1]
        if lab == "WT"
            data = readdlm(file, ',', Int)[1:end,1]
            first_idx = argmax(data)
            end_idx = min(first_idx + 45, length(data))
        elseif lab == "rbmA"
            data = readdlm(file, ',', Int)[1:end,1] .+ 0.01
            data_max = maximum(data)
            fc = data ./ data_max
            min_fc = minimum(fc)
            @show min_fc
            fc = fc .- min_fc
            end_idx = findfirst(x->x<0.001, fc) 
            first_idx = end_idx - 45 
            @show fc[first_idx:end_idx]
        end
        data_subset = data[first_idx:end_idx] 
        data = data_subset
        data_norm =  data ./ data[1]
        data = data_norm
        if lab == "WT"
            color = 1 
            colormap = :Pastel1_5
            if WT_seen
                condition = nothing 
            else
                WT_seen = true
                condition = "Wild-type"
            end
        elseif lab == "rbmA"
            color =5 
            colormap = :Pastel1_5
            if rbmA_seen
                condition = nothing 
            else
                rbmA_seen = true
                condition = rich("Δ", rich("rbmA"; font=:italic))
            end
        end
        xs = 0:6:length(data)-1
        xaxis = 1:length(data)
        xaxis /= 6
        lines!(ax, xaxis, data, color=color, colormap=colormap, colorrange = (1, 5), label=condition)
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
    save(plots_folder*"/WT_rbmA_bulk.svg", fig)
end

main()
