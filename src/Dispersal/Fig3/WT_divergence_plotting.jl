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
                      "/mnt/h/Dispersal/WT_replicate5_processed/Displacements/"]
    dx = 8
    dy = 8
    dz = 8
    WT_seen = false
    cheY_seen = false
    rbmB_seen = false
    lapG_seen = false
    bulk_files = [f for f in readdir(plots_folder, join=true) if occursin("processed.csv", f)]
    WT_dispersal = []
    WT_growth = []
    conditions = []
    for vector_folder in vector_folders
        bulk_file = [f for f in bulk_files if occursin(split(vector_folder, "/")[end-2], f)][1]
        bulk_data = readdlm(bulk_file, ',', Int)[1:end,1]
        first_idx = argmax(bulk_data)
        end_idx = min(first_idx+45, length(bulk_data))
        piv_files = [f for f in readdir(vector_folder, join=true) if occursin("piv", f)]
        dispersal_divergences = Array{Float64, 1}(undef, end_idx-first_idx+1)
        growth_divergences = Array{Float64, 1}(undef, first_idx-1)
        for i in 1:first_idx-1 
            u, v, w, flags = load(piv_files[i], "u", "v", "w", "flags")
            ui = permutedims(u, [2,1,3])
            vi = permutedims(v, [2,1,3])
            wi = permutedims(w, [2,1,3])
            flags = permutedims(flags, [2,1,3])
            ui[flags .> 0] .= NaN
            vi[flags .> 0] .= NaN
            wi[flags .> 0] .= NaN
            net_divergence = calculate_divergence(ui, vi, wi, dx, dy, dz)
            growth_divergences[i] = net_divergence
        end
        growth_divergences = mean(growth_divergences) 
        push!(WT_growth, growth_divergences)
        push!(conditions, "Growth")
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
            dispersal_divergences[i-first_idx+1] = net_divergence
        end
        dispersal_divergences = mean(dispersal_divergences) 
        push!(WT_dispersal, dispersal_divergences)
        push!(conditions, "Dispersal")
    end
    data = vcat(WT_growth, WT_dispersal)
    category_num = Int.(repeat(1:2, inner = 5))
    fig = Figure(size=(2.5*72, 3*72))
    ax = Axis(fig[1, 1])
    boxplot!(ax, category_num, Float64.(data); show_outliers=false, color="#fbb4ae")
    ax.xticks=(1:2, unique(conditions))
    ax.xticklabelrotation=45
    ax.xlabel=""
    ax.ylabel=rich("Divergence rate (h", superscript("-1"), ")")
    ax.title=""
    ax.rightspinevisible = false
    ax.topspinevisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    save(plots_folder*"/WT_convergence.svg", fig)
end

main()
