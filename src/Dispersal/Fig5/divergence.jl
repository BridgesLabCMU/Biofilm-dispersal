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
using NaNStatistics
using ImageMorphology
using HypothesisTests

function choose_label(folder)
    if occursin("WT", folder)
        return "WT"
    elseif occursin("rbmA", folder)
        return "rbmA"
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
    conditions = ["Wild-type", L"$\Delta rbmA$"]
    images_folder = "/mnt/h/Dispersal/WT_replicate1_processed/"
    plots_folder = "/mnt/h/Dispersal/Plots"
    mask_files = sort([f for f in readdir(images_folder, join=true) if occursin("mask_isotropic", f)], lt=natural)
    mask_t1 = load(mask_files[1])
    vector_folders = ["/mnt/h/Dispersal/WT_replicate1_processed/Displacements/", 
                      "/mnt/h/Dispersal/WT_replicate2_processed/Displacements/", 
                      "/mnt/h/Dispersal/WT_replicate3_processed/Displacements/", 
                      "/mnt/h/Dispersal/WT_replicate4_processed/Displacements/", 
                      "/mnt/h/Dispersal/WT_replicate5_processed/Displacements/", 
                      "/mnt/h/Dispersal/rbmA_replicate1_processed/Displacements/", 
                      "/mnt/h/Dispersal/rbmA_replicate2_processed/Displacements/", 
                      "/mnt/h/Dispersal/rbmA_replicate3_processed/Displacements/", 
                      "/mnt/h/Dispersal/rbmA_replicate4_processed/Displacements/", 
                      "/mnt/h/Dispersal/rbmA_replicate5_processed/Displacements/"]
    images_folders = ["/mnt/h/Dispersal/WT_replicate1_processed/", 
                      "/mnt/h/Dispersal/WT_replicate2_processed/", 
                      "/mnt/h/Dispersal/WT_replicate3_processed/", 
                      "/mnt/h/Dispersal/WT_replicate4_processed/", 
                      "/mnt/h/Dispersal/WT_replicate5_processed/", 
                      "/mnt/h/Dispersal/rbmA_replicate1_processed/", 
                      "/mnt/h/Dispersal/rbmA_replicate2_processed/", 
                      "/mnt/h/Dispersal/rbmA_replicate3_processed/", 
                      "/mnt/h/Dispersal/rbmA_replicate4_processed/", 
                      "/mnt/h/Dispersal/rbmA_replicate5_processed/"]
    dx = 8
    dy = 8
    dz = 8
    WT_seen = false
    rbmA_seen = false
    bulk_files = [f for f in readdir(plots_folder, join=true) if occursin("processed.csv", f)]
    logocolors = Colors.JULIA_LOGO_COLORS
    WT_averages = []
    rbmA_averages = []
    conditions = []
    for (j, vector_folder) in enumerate(vector_folders)
        mask_files = sort([f for f in readdir(images_folders[j], join=true) if occursin("mask_isotropic", f)], lt=natural)
        mask_t1 = load(mask_files[1])
        bulk_file = [f for f in bulk_files if occursin(split(vector_folder, "/")[end-2], f)][1]
        if occursin("WT", vector_folder)
            bulk_data = readdlm(bulk_file, ',', Int)[1:end,1]
            first_idx = argmax(bulk_data)
            end_idx = min(first_idx+45, length(bulk_data))
        elseif occursin("rbmA", vector_folder)
            bulk_data = readdlm(bulk_file, ',', Int)[1:end,1] .+ 0.01
            data_max = maximum(bulk_data)
            fc = bulk_data ./ data_max
            min_fc = minimum(fc)
            fc = fc .- min_fc
            end_idx = findfirst(x->x<0.001, fc) 
            first_idx = 6 #end_idx - 45 
        end
        piv_files = sort([f for f in readdir(vector_folder, join=true) if occursin("piv", f)], lt=natural)
        if end_idx > length(piv_files)
            end_idx = length(piv_files)  
        end
        dummy_u = load(piv_files[first_idx], "u")
        center_of_mass = compute_com(mask_t1, dummy_u)
        height, width, depth = size(dummy_u)
        u_tot = zeros(Float32, width, height, depth)
        v_tot = zeros(Float32, width, height, depth)
        w_tot = zeros(Float32, width, height, depth) 
        for i in first_idx:end_idx 
            u, v, w, flags = load(piv_files[i], "u", "v", "w", "flags")
            u[flags .> 0] .= 0
            v[flags .> 0] .= 0
            w[flags .> 0] .= 0
            ui = permutedims(u, [2,1,3])
            vi = permutedims(v, [2,1,3])
            wi = permutedims(w, [2,1,3])
            u_tot = u_tot + ui
            v_tot = v_tot + vi
            w_tot = w_tot + wi
        end
        divergences = calculate_radial_component(u_tot, v_tot, w_tot.*4, dx, dy, dz, center_of_mass)
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
        elseif lab == "rbmA"
            push!(rbmA_averages, divergences)
            c = logocolors.green
            if rbmA_seen
                condition = ""
            else
                condition = latexstring("\\Delta $(lab)")
                rbmA_seen = true
                push!(conditions, rich("Δ", rich("rbmA"; font=:italic)))
            end
		end
    end
    data = vcat(WT_averages, rbmA_averages) .* 0.065
    writedlm("$(plots_folder)/Fig5D.csv", data, ",")                                                                 
    @show pvalue(UnequalVarianceTTest(Float64.(WT_averages), Float64.(rbmA_averages)))
    averages = [mean(data[1:5]), mean(data[6:10])]
    maxes = [maximum(data[1:5]), maximum(data[6:10])]
    mins = [minimum(data[1:5]), minimum(data[6:10])]
    category_num_swarm = Int.(repeat(1:2, inner = 5))
    fig = Figure(size=(3*72, 3.5*72))
	category_num = Int.(1:2)
	category_num_swarm = Int.(repeat(1:2, inner=5))
	ax = Axis(fig[1, 1])
    colormap1 = [[:black]; Makie.wong_colors()[3]]
    colormap2 = [[:white]; Makie.wong_colors()[3]]
	crossbar!(ax, category_num, averages, mins, maxes; 
			  color=:white, midlinecolor=colormap1, colormap1, colorrange=(1,2))
    plt = beeswarm!(ax, category_num_swarm, Float64.(data), color = category_num_swarm, algorithm=UniformJitter(), strokewidth=1)
	plt.colormap[] = colormap2 
    ax.xticks=(1:2, conditions)
    ax.xticklabelrotation=45
    ax.xlabel=""
    ax.ylabel="Radial displacement (µm)"
    ax.title=""
    ax.rightspinevisible = false
    ax.topspinevisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    save(plots_folder*"/WT_rbmA_convergence.svg", fig)
end
main()
