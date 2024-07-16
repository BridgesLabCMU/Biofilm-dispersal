using FileIO 
using Makie 
using GLMakie
using CairoMakie
using SwarmMakie
CairoMakie.activate!(type = "pdf")
set_theme!(fonts = (; regular = "TeX Gyre Heros Regular", bold = "TeX Gyre Heros Regular"))
using NaturalSort
using StatsBase
using NaNStatistics
using ImageMorphology
using Images
using HistogramThresholding
using DelimitedFiles
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
    return radial_components
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
	return center_of_mass_vectors, center_of_mass_images
end

function main()
    images_folders = ["/mnt/h/Dispersal/WT_replicate1_processed", "/mnt/h/Dispersal/WT_replicate2_processed", "/mnt/h/Dispersal/WT_replicate3_processed", "/mnt/h/Dispersal/WT_replicate4_processed", "/mnt/h/Dispersal/WT_replicate5_processed",
                     "/mnt/h/Dispersal/rbmA_replicate1_processed", "/mnt/h/Dispersal/rbmA_replicate2_processed", "/mnt/h/Dispersal/rbmA_replicate3_processed", "/mnt/h/Dispersal/rbmA_replicate4_processed", "/mnt/h/Dispersal/rbmA_replicate5_processed"]
    plots_folder = "/mnt/h/Dispersal/Plots"
    WT_averages = []
    rbmA_averages = []
    conditions = []
    WT_seen = false
    rbmA_seen = false
    for images_folder in images_folders
        folder = images_folder*"/Displacements/"
        files = sort([f for f in readdir(folder, join=true) if occursin("piv", f)], lt=natural)
        mask_files = sort([f for f in readdir(images_folder, join=true) if occursin("mask_isotropic", f)], lt=natural)
        mask_t1 = load(mask_files[1])
        net = readdlm(plots_folder*"/"*basename(images_folder)*".csv", ',', Int)[1:end,1]
        first_index = argmax(net)
        end_index = first_index + 20
        first_index = max(1, first_index - 4)

        x, y, z, u_dummy = load(files[1], "x", "y", "z", "u")
        x = Int.(x)
        y = Int.(y)
        z = Int.(z)

        flags_tot = zeros(Bool, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])
        u_plot = zeros(Bool, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])
        v_plot = zeros(Bool, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])
        w_plot = zeros(Bool, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])

        n_out = Array{Float32, 1}(undef, length(files))
        slice=round(Int, 0.2*size(u_dummy)[3])

        # Calculate cumulative displacements
        for i in first_index:end_index
            u, v, w, flags = load(files[i], "u", "v", "w", "flags")
            u[flags .> 0] .= 0
            v[flags .> 0] .= 0
            w[flags .> 0] .= 0
            u = mapwindow(median, u, (3,3,3))
            v = mapwindow(median, v, (3,3,3))
            w = mapwindow(median, w, (3,3,3))
            u_plot += permutedims(u, (2,1,3))
            v_plot += permutedims(v, (2,1,3))
            w_plot += permutedims(w, (2,1,3))
        end
    
        center_of_mass, com_images = compute_com(mask_t1, u_plot)
        radial_component = calculate_radial_component(u_plot, v_plot, w_plot.*4, 8, 8, 8, center_of_mass)

        # Calculate local mean, std
        u_mean = mapwindow(mean, u_plot, (9,9,9))
        v_mean = mapwindow(mean, v_plot, (9,9,9))
        w_mean = mapwindow(mean, w_plot, (9,9,9))
        u_std = mapwindow(std, u_plot, (9,9,9))
        v_std = mapwindow(std, v_plot, (9,9,9))
        w_std = mapwindow(std, w_plot, (9,9,9))
        radial_component_mean = calculate_radial_component(u_mean, v_mean, w_mean.*4, 8, 8, 8, center_of_mass)
        radial_component_std = calculate_radial_component(u_std, v_std, w_std.*4, 8, 8, 8, center_of_mass)
        
        # Calculate distance
        distance = sqrt.(u_plot.^2 .+ v_plot.^2 .+ w_plot.^2)
        mean_distance = sqrt.(u_mean.^2 .+ v_mean.^2 .+ w_mean.^2)

        channel_volume = sum(radial_component .> radial_component_mean .+ 4) / sum(distance .>= 0)
        lab = choose_label(folder)
        if lab == "WT"
            push!(WT_averages, channel_volume)
            if WT_seen
                condition = ""
            else
                condition = "Wild-type"
                WT_seen = true
                push!(conditions, "Wild-type")
            end
        elseif lab == "rbmA"
            push!(rbmA_averages, channel_volume)
            c = :coral2 
            if rbmA_seen
                condition = ""
            else
                rbmA_seen = true
                push!(conditions, rich("Δ", rich("rbmA"; font=:italic)))
            end
        end
    end
    data = vcat(WT_averages, rbmA_averages)
    writedlm("$(plots_folder)/Fig5E.csv", data, ",")
    @show pvalue(UnequalVarianceTTest(Float64.(WT_averages), Float64.(rbmA_averages)))
    averages = [mean(data[1:5]), mean(data[6:10])]
    maxes = [maximum(data[1:5]), maximum(data[6:10])]
    mins = [minimum(data[1:5]), minimum(data[6:10])]
    category_num_swarm = Int.(repeat(1:2, inner = 5))
    fig = Figure(size=(3*72, 3.5*72))
	category_num = Int.(1:2)
	category_num_swarm = Int.(repeat(1:2, inner=5))
	ax = CairoMakie.Axis(fig[1, 1])
    colormap1 = [[:black]; Makie.wong_colors()[4]]
    colormap2 = [[:white]; Makie.wong_colors()[4]]
	crossbar!(ax, category_num, averages, mins, maxes; 
			  color=:white, midlinecolor=colormap1, colormap1, colorrange=(1,2))
    plt = beeswarm!(ax, category_num_swarm, Float64.(data), color = category_num_swarm, algorithm=UniformJitter(), strokewidth=1)
	plt.colormap[] = colormap2 
    ax.xticks=(1:2, conditions)
    ax.xticklabelrotation=45
    ax.xlabel=""
    ax.ylabel=rich("Channel volume fraction (a.u.)")
    ax.title=""
    ax.rightspinevisible = false
    ax.topspinevisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    ylims!(ax, 0.0, nothing)
    save(plots_folder*"/channels_rbmA.svg", fig)
end
main()
