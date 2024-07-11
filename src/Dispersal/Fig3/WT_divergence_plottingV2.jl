using Makie 
using GLMakie
using CairoMakie
CairoMakie.activate!(type = "svg")
using SwarmMakie
using LaTeXStrings
using StatsBase
using DelimitedFiles
using Colors
using NaturalSort
using FileIO
using Makie.Colors
using ImageMorphology
using NaNStatistics
using HypothesisTests

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

function calculate_radial_component(u, v, w, dx,dy,dz, central_point)
    nx, ny, nz = size(u)
    x = reshape(1:nx, nx, 1, 1)
    y = reshape(1:ny, 1, ny, 1)
    z = reshape(1:nz, 1, 1, nz)
    dx = x .- central_point[1]
    dy = y .- central_point[2]
    dz = z .- central_point[3]
	dx = repeat(dx, 1, ny, nz)
    dy = repeat(dy, nx, 1, nz)
    dz = repeat(dz, nx, ny, 1)
    magnitudes = sqrt.(dx.^2 + dy.^2 + dz.^2)
    dx_norm = dx ./ magnitudes
    dy_norm = dy ./ magnitudes
    dz_norm = dz ./ magnitudes
    radial_components = u .* dx_norm + v .* dy_norm + w .* dz_norm
    return nanmean(radial_components)
end

function compute_com(mask, vector_field)
	mask = label_components(mask)
	areas = component_lengths(mask)
	max_label = argmax(areas[1:end])
	mask[mask .!= max_label] .= 0
	mask[mask .== max_label] .= 1 
    center_of_mass_images = component_centroids(mask)[1]
    relative_com = [x for x in center_of_mass_images] ./ [x for x in size(mask)] 
    center_of_mass_images = [center_of_mass_images[2], center_of_mass_images[1], center_of_mass_images[3]]
    relative_com = [relative_com[2], relative_com[1], relative_com[3]]
    center_of_mass_vectors = relative_com .* [x for x in size(vector_field)] .- 1
    center_of_mass_vectors[3] = 0
	return center_of_mass_vectors
end

function main()
    conditions = ["Wild-type", L"$\Delta cheY$",L"$\Delta lapG$",L"$\Delta rbmB$"]
    images_folders = ["/mnt/h/Dispersal/WT_replicate1_processed/", 
                      "/mnt/h/Dispersal/WT_replicate2_processed/", 
                      "/mnt/h/Dispersal/WT_replicate3_processed/", 
                      "/mnt/h/Dispersal/WT_replicate4_processed/", 
                      "/mnt/h/Dispersal/WT_replicate5_processed/"]
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
    for (j, vector_folder) in enumerate(vector_folders)
        mask_files = sort([f for f in readdir(images_folders[j], join=true) if occursin("mask_isotropic", f)], lt=natural)
        mask_t1 = load(mask_files[1])
        bulk_file = [f for f in bulk_files if occursin(split(vector_folder, "/")[end-2], f)][1]
        bulk_data = readdlm(bulk_file, ',', Int)[1:end,1]
        first_idx = argmax(bulk_data)
        end_idx = min(first_idx+45, length(bulk_data))
        piv_files = sort([f for f in readdir(vector_folder, join=true) if occursin("piv", f)], lt=natural)
        dispersal_divergences = Array{Float64, 1}(undef, end_idx-first_idx+1)
        growth_divergences = Array{Float64, 1}(undef, first_idx-1)
        dummy_u = load(piv_files[first_idx], "u")
        height, width, depth = size(dummy_u)
        center_of_mass = compute_com(mask_t1, dummy_u)
        u_dispersal = zeros(Float32, width, height, depth)
        v_dispersal = zeros(Float32, width, height, depth)
        w_dispersal = zeros(Float32, width, height, depth) 
        u_growth = zeros(Float32, width, height, depth)
        v_growth = zeros(Float32, width, height, depth)
        w_growth = zeros(Float32, width, height, depth)
        for i in 1:first_idx-1 
            u, v, w, flags = load(piv_files[i], "u", "v", "w", "flags")
            u[flags .> 0] .= 0
            v[flags .> 0] .= 0
            w[flags .> 0] .= 0
            ui = permutedims(u, [2,1,3])
            vi = permutedims(v, [2,1,3])
            wi = permutedims(w, [2,1,3])
            u_growth = u_growth + ui
            v_growth = v_growth + vi
            w_growth = w_growth + wi
        end
        growth_divergences = calculate_radial_component(u_growth, v_growth, w_growth.*4, 8, 8, 8, center_of_mass) * 0.065 
        push!(WT_growth, growth_divergences)
        push!(conditions, "Growth")
        for i in first_idx:end_idx 
            u, v, w, flags = load(piv_files[i], "u", "v", "w", "flags")
            u[flags .> 0] .= 0
            v[flags .> 0] .= 0
            w[flags .> 0] .= 0
            ui = permutedims(u, [2,1,3])
            vi = permutedims(v, [2,1,3])
            wi = permutedims(w, [2,1,3])
            u_dispersal = u_dispersal + ui
            v_dispersal = v_dispersal + vi
            w_dispersal = w_dispersal + wi
        end
        dispersal_divergences = calculate_radial_component(u_dispersal, v_dispersal, w_dispersal.*4, 8, 8, 8, center_of_mass) * 0.065
        push!(WT_dispersal, dispersal_divergences)
        push!(conditions, "Dispersal")
    end
    data = vcat(WT_growth, WT_dispersal)
    @show pvalue(UnequalVarianceTTest(Float64.(WT_growth), Float64.(WT_dispersal)))
    averages = [mean(data[1:5]), mean(data[6:10])]
    maxes = [maximum(data[1:5]), maximum(data[6:10])]
    mins = [minimum(data[1:5]), minimum(data[6:10])]
    category_num = Int.(1:2)
    category_num_swarm = Int.(repeat(1:2, inner = 5))
    fig = Figure(size=(3*72, 3*72))
    ax = CairoMakie.Axis(fig[1, 1])
	crossbar!(ax, category_num, averages, mins, maxes; 
			  color=:white, midlinecolor=:black)
    plt = beeswarm!(ax, category_num_swarm, Float64.(data), strokecolor=:black, color=:white, algorithm=UniformJitter(), strokewidth=1)
    ax.xticks=(1:2, unique(conditions))
    ax.xticklabelrotation=45
    ax.xlabel=""
    ax.ylabel="Radial displacement (µm)"
    ax.title=""
    ax.rightspinevisible = false
    ax.topspinevisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    save(plots_folder*"/WT_convergence.svg", fig)
end

main()
