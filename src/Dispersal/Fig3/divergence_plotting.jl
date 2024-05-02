using Makie 
using GLMakie
using CairoMakie
CairoMakie.activate!(type = "svg")
using LaTeXStrings
using StatsBase
using DelimitedFiles
using Colors
using NaturalSort
using FileIO
using Makie.Colors

function choose_label(folder)
    if occursin("WT", folder)
        return "WT"
    elseif occursin("cheY", folder)
        return "cheY"
    elseif occursin("lapG", folder)
        return "lapG"
    elseif occursin("rbmB", folder)
        return "rbmB"
    end
end

function calculate_divergence(u, v, w, dx, dy, dz)
    nx, ny, nz = size(u)
    net_divergence = 0.0
    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                if !isnan(u[i,j,k])
                    dudx = (i==nx && !isnan(u[i-1,j,k])) ? (u[i,j,k]-u[i-1,j,k])/dx : (i==1 && !isnan(u[i+1,j,k])) ? (u[i+1,j,k]-u[i,j,k])/dx : (i!=nx && i!=1 && !isnan(u[i-1,j,k]) && !isnan(u[i+1,j,k])) ? (u[i+1,j,k]-u[i-1,j,k])/(2dx) : NaN 
				else
                    dudx = NaN 
                end
                if !isnan(v[i,j,k])
                    dvdy = (j==ny && !isnan(v[i,j-1,k])) ? (v[i,j,k]-v[i,j-1,k])/dx : (j==1 && !isnan(v[i,j+1,k])) ? (v[i,j+1,k]-v[i,j,k])/dx : (j!=ny && j!=1 && !isnan(v[i,j-1,k]) && !isnan(v[i,j+1,k])) ? (v[i,j+1,k]-v[i,j-1,k])/(2dx) : NaN 
                else
                    dvdy = NaN 
                end
                if !isnan(w[i,j,k])
                    dwdz = (k==nz && !isnan(w[i,j,k-1])) ? (w[i,j,k]-w[i,j,k-1])/dx : (k==1 && !isnan(w[i,j,k+1])) ? (w[i,j,k+1]-w[i,j,k])/dx : (k!=nz && k!=1 && !isnan(w[i,j,k-1]) && !isnan(w[i,j,k+1])) ? (w[i,j,k+1]-w[i,j,k-1])/(2dx) : NaN 
                else
                    dwdz = NaN 
                end
                if !isnan(dudx + dvdy + dwdz)
                    net_divergence += dudx + dvdy + dwdz
                end
            end
        end
    end
    return net_divergence
end

function main()
    conditions = ["Wild-type", L"$\Delta cheY$",L"$\Delta lapG$",L"$\Delta rbmB$"]
    plots_folder = "/mnt/h/Dispersal/Plots"
    vector_folders = ["/mnt/h/Dispersal/WT_replicate1_processed/Displacements/", 
                      "/mnt/h/Dispersal/WT_replicate2_processed/Displacements/", 
                      "/mnt/h/Dispersal/WT_replicate3_processed/Displacements/", 
                      "/mnt/h/Dispersal/WT_replicate4_processed/Displacements/", 
                      "/mnt/h/Dispersal/WT_replicate5_processed/Displacements/", 
                      "/mnt/h/Dispersal/cheY_replicate1_processed/Displacements/", 
                      "/mnt/h/Dispersal/cheY_replicate2_processed/Displacements/", 
                      "/mnt/h/Dispersal/cheY_replicate3_processed/Displacements/", 
                      "/mnt/h/Dispersal/cheY_replicate4_processed/Displacements/", 
                      "/mnt/h/Dispersal/cheY_replicate5_processed/Displacements/", 
                      "/mnt/h/Dispersal/lapG_replicate1_processed/Displacements/", 
                      "/mnt/h/Dispersal/lapG_replicate2_processed/Displacements/", 
                      "/mnt/h/Dispersal/lapG_replicate3_processed/Displacements/", 
                      "/mnt/h/Dispersal/lapG_replicate4_processed/Displacements/", 
                      "/mnt/h/Dispersal/lapG_replicate5_processed/Displacements/", 
                      "/mnt/h/Dispersal/rbmB_replicate1_processed/Displacements/", 
                      "/mnt/h/Dispersal/rbmB_replicate2_processed/Displacements/", 
                      "/mnt/h/Dispersal/rbmB_replicate3_processed/Displacements/", 
                      "/mnt/h/Dispersal/rbmB_replicate4_processed/Displacements/", 
                      "/mnt/h/Dispersal/rbmB_replicate5_processed/Displacements/"]
    dx = 8
    dy = 8
    dz = 8
    WT_seen = false
    cheY_seen = false
    rbmB_seen = false
    lapG_seen = false
    bulk_files = [f for f in readdir(plots_folder, join=true) if occursin("processed.csv", f)]
    logocolors = Colors.JULIA_LOGO_COLORS
    WT_averages = []
    cheY_averages = []
    rbmB_averages = []
    lapG_averages = []
    conditions = []
    for vector_folder in vector_folders
        bulk_file = [f for f in bulk_files if occursin(split(vector_folder, "/")[end-2], f)][1]
        bulk_data = readdlm(bulk_file, ',', Int)[1:end,1]
        first_idx = argmax(bulk_data)
        end_idx = min(first_idx+45, length(bulk_data))
        piv_files = [f for f in readdir(vector_folder, join=true) if occursin("piv", f)]
        divergences = Array{Float64, 1}(undef, end_idx-first_idx+1)
        for i in first_idx:end_idx 
            u, v, w, flags = load(piv_files[i], "u", "v", "w", "flags")
            ui = permutedims(u, [2,1,3])
            vi = permutedims(v, [2,1,3])
            wi = permutedims(w, [2,1,3])
            flags = permutedims(flags, [2,1,3])
            ui[flags .> 0] .= NaN
            vi[flags .> 0] .= NaN
            wi[flags .> 0] .= NaN
            net_divergence = calculate_divergence(ui, vi, wi, dx, dy, dz)
            divergences[i-first_idx+1] = net_divergence
        end
        divergences = mean(divergences) 
        lab = choose_label(vector_folder)
        if lab == "WT"
            push!(WT_averages, divergences)
            c = logocolors.blue
            if WT_seen
                condition = ""
            else
                condition = "Wild-type"
                WT_seen = true
                push!(conditions, "Wild-type")
            end
        elseif lab == "cheY"
            push!(cheY_averages, divergences)
            c = logocolors.green
            if cheY_seen
                condition = ""
            else
                condition = latexstring("\\Delta $(lab)")
                cheY_seen = true
                push!(conditions, rich("Δ", rich("cheY"; font=:italic)))
            end
        elseif lab == "rbmB"
            push!(rbmB_averages, divergences)
            c = :coral2 
            if rbmB_seen
                condition = ""
            else
                condition = latexstring("\\Delta $(lab)")
                rbmB_seen = true
                push!(conditions, rich("Δ", rich("rbmB"; font=:italic)))
            end
        elseif lab == "lapG"
            push!(lapG_averages, divergences)
            c = logocolors.purple
            if lapG_seen
                condition = ""
            else
                condition = latexstring("\\Delta $(lab)")
                lapG_seen = true
                push!(conditions, rich("Δ", rich("lapG"; font=:italic)))
            end
        end
    end
    data = vcat(WT_averages, cheY_averages, lapG_averages, rbmB_averages)
    category_num = Int.(repeat(1:4, inner = 5))
    fig = Figure(size=(3*72, 3*72))
    ax = Axis(fig[1, 1])
    colormap = Makie.to_colormap(:Pastel1_4)
    boxplot!(ax, category_num, Float64.(data); show_outliers=false, color=map(category_num->mod1(category_num,4),category_num), colormap, colorrange=(1,4))
    ax.xticks=(1:4, conditions)
    ax.xticklabelrotation=45
    ax.xlabel=""
    ax.ylabel=rich("Divergence rate (h", superscript("-1"), ")")
    ax.title=""
    ax.rightspinevisible = false
    ax.topspinevisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    save(plots_folder*"/all_convergence.svg", fig)
end

main()
