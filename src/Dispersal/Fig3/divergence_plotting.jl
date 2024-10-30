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
    elseif occursin("cheY", folder)
        return "cheY"
    elseif occursin("lapG", folder)
        return "lapG"
    elseif occursin("rbmB", folder)
        return "rbmB"
    end
end

function calculate_radial_component(u, v, w, dx,dy,dz, central_point)
    # Vector field comes in in units of voxels
    # Function to calculate radial components, given PIV vector field, 
    nx, ny, nz = size(u)
    # Get grid coordinates (in grid units)
    x = reshape(1:nx, nx, 1, 1)
    y = reshape(1:ny, 1, ny, 1)
    z = reshape(1:nz, 1, 1, nz)
    # Grid coordinates relative to the center of mass
    dx = x .- central_point[1]
    dy = y .- central_point[2]
    dz = z .- central_point[3]
	dx = repeat(dx, 1, ny, nz)
    dy = repeat(dy, nx, 1, nz)
    dz = repeat(dz, nx, ny, 1)
    # Magnitude of the displacement
    magnitudes = sqrt.(dx.^2 + dy.^2 + dz.^2)
    dx_norm = dx ./ magnitudes
    dy_norm = dy ./ magnitudes
    dz_norm = dz ./ magnitudes
    # Radial component (units of voxels)
    radial_components = u .* dx_norm + v .* dy_norm + w .* dz_norm
    return radial_components
end

function compute_com(mask, vector_field)
	mask = label_components(mask)
	areas = component_lengths(mask)
    # Biofilm = largest connected component
	max_label = argmax(areas[1:end])
	mask[mask .!= max_label] .= 0
	mask[mask .== max_label] .= 1 
    center_of_mass_images = component_centroids(mask)[1]
    relative_com = [x for x in center_of_mass_images] ./ [x for x in size(mask)] 
    # Switch to physical coordinates
    center_of_mass_images = [center_of_mass_images[2], center_of_mass_images[1], center_of_mass_images[3]]
    relative_com = [relative_com[2], relative_com[1], relative_com[3]]
    # Convert center of mass to PIV coordinates
    center_of_mass_vectors = relative_com .* [x for x in size(vector_field)] .- 1
    # Set z = 0 (center of mass at the bottom of the biofilm)
    center_of_mass_vectors[3] = 0
	return center_of_mass_vectors, center_of_mass_images
end

function main()
    conditions = ["Wild-type", L"$\Delta cheY$",L"$\Delta lapG$",L"$\Delta rbmB$"]
    images_folder = "/mnt/h/Dispersal/WT_replicate1_processed/"
    plots_folder = "/mnt/h/Dispersal/Plots"
    mask_files = sort([f for f in readdir(images_folder, join=true) if occursin("mask_isotropic", f)], lt=natural)
    mask_t1 = load(mask_files[1])
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
    images_folders = ["/mnt/h/Dispersal/WT_replicate1_processed/", 
                      "/mnt/h/Dispersal/WT_replicate2_processed/", 
                      "/mnt/h/Dispersal/WT_replicate3_processed/", 
                      "/mnt/h/Dispersal/WT_replicate4_processed/", 
                      "/mnt/h/Dispersal/WT_replicate5_processed/", 
                      "/mnt/h/Dispersal/cheY_replicate1_processed/", 
                      "/mnt/h/Dispersal/cheY_replicate2_processed/", 
                      "/mnt/h/Dispersal/cheY_replicate3_processed/", 
                      "/mnt/h/Dispersal/cheY_replicate4_processed/", 
                      "/mnt/h/Dispersal/cheY_replicate5_processed/", 
                      "/mnt/h/Dispersal/rbmB_replicate1_processed/", 
                      "/mnt/h/Dispersal/rbmB_replicate2_processed/", 
                      "/mnt/h/Dispersal/rbmB_replicate3_processed/", 
                      "/mnt/h/Dispersal/rbmB_replicate4_processed/", 
                      "/mnt/h/Dispersal/rbmB_replicate5_processed/",
                      "/mnt/h/Dispersal/lapG_replicate1_processed/", 
                      "/mnt/h/Dispersal/lapG_replicate2_processed/", 
                      "/mnt/h/Dispersal/lapG_replicate3_processed/", 
                      "/mnt/h/Dispersal/lapG_replicate4_processed/", 
                      "/mnt/h/Dispersal/lapG_replicate5_processed/"] 
    dx = 8
    dy = 8
    dz = 8
    WT_seen = false
    cheY_seen = false
    rbmB_seen = false
    lapG_seen = false
    bulk_files = [f for f in readdir(plots_folder, join=true) if occursin("processed.csv", f)]
    WT_averages = []
    cheY_averages = []
    rbmB_averages = []
    lapG_averages = []
    conditions = []
    for (j, vector_folder) in enumerate(vector_folders)
        mask_files = sort([f for f in readdir(images_folders[j], join=true) if occursin("mask_isotropic", f)], lt=natural)
        mask_t1 = load(mask_files[1])
        bulk_file = [f for f in bulk_files if occursin(split(vector_folder, "/")[end-2], f)][1]
        bulk_data = readdlm(bulk_file, ',', Int)[1:end,1]
        first_idx = argmax(bulk_data)
        end_idx = min(first_idx+45, length(bulk_data))
        piv_files = sort([f for f in readdir(vector_folder, join=true) if occursin("piv", f)], lt=natural)
        dummy_u = load(piv_files[first_idx], "u")
        center_of_mass = compute_com(mask_t1, dummy_u)
        height, width, depth = size(dummy_u)
        u_tot = zeros(Float32, width, height, depth)
        v_tot = zeros(Float32, width, height, depth)
        w_tot = zeros(Float32, width, height, depth) 
        # Calculate total displacement
        for i in first_idx:end_idx 
            u, v, w, flags = load(piv_files[i], "u", "v", "w", "flags")
            # Set spurious vectors to 0
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
            if WT_seen
                condition = ""
            else
                condition = "Wild-type"
                WT_seen = true
                push!(conditions, "Wild-type")
            end
        elseif lab == "cheY"
            push!(cheY_averages, divergences)
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
            if lapG_seen
                condition = ""
            else
                condition = latexstring("\\Delta $(lab)")
                lapG_seen = true
                push!(conditions, rich("Δ", rich("lapG"; font=:italic)))
            end
        end
    end
    data = vcat(WT_averages, cheY_averages, rbmB_averages, lapG_averages) .* 0.065 # Convert from image voxels to µm
    writedlm("$(plots_folder)/Fig4F.csv", data, ",")
    @show pvalue(UnequalVarianceTTest(Float64.(WT_averages), Float64.(cheY_averages)))
    @show pvalue(UnequalVarianceTTest(Float64.(WT_averages), Float64.(lapG_averages)))
    @show pvalue(UnequalVarianceTTest(Float64.(WT_averages), Float64.(rbmB_averages)))
    averages = [mean(data[1:5]), mean(data[6:10]), mean(data[11:15]), mean(data[16:20])]
    maxes = [maximum(data[1:5]), maximum(data[6:10]), maximum(data[11:15]), maximum(data[16:20])]
    mins = [minimum(data[1:5]), minimum(data[6:10]), minimum(data[11:15]), minimum(data[16:20])]
    category_num_swarm = Int.(repeat(1:4, inner = 5))
    fig = Figure(size=(4*72, 3.5*72))
	category_num = Int.(1:4)
	category_num_swarm = Int.(repeat(1:4, inner=5))
	ax = Axis(fig[1, 1])
    colormap1 = [[:black]; Makie.wong_colors()[1:2]; Makie.wong_colors()[4]]
    colormap2 = [[:white]; Makie.wong_colors()[1:2]; Makie.wong_colors()[4]]
	crossbar!(ax, category_num, averages, mins, maxes; 
			  color=:white, midlinecolor=colormap1, colormap1, colorrange=(1,4))
    plt = beeswarm!(ax, category_num_swarm, Float64.(data), color = category_num_swarm, algorithm=UniformJitter(), strokewidth=1)
    plt.colormap[] = colormap2 
    ax.xticks=(1:4, conditions)
    ax.xticklabelrotation=45
    ax.xlabel=""
    ax.ylabel="Radial displacement (µm)"
    ax.title=""
    ax.rightspinevisible = false
    ax.topspinevisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    save(plots_folder*"/all_convergence.svg", fig)
end

main()
