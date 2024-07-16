using FileIO 
using Images
using Makie 
using GLMakie
using CairoMakie
CairoMakie.activate!(type = "pdf")
set_theme!(fonts = (; regular = "TeX Gyre Heros Regular", bold = "TeX Gyre Heros Regular"))
using NaturalSort
using StatsBase
using NaNStatistics
using ImageMorphology

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
    folder = "/mnt/h/Dispersal/WT_replicate1_processed/Displacements/"
    images_folder = "/mnt/h/Dispersal/WT_replicate1_processed/"
    plots_folder = "/mnt/h/Dispersal/Plots"
    mask_files = sort([f for f in readdir(images_folder, join=true) if occursin("mask_isotropic", f)], lt=natural)
    mask_t1 = load(mask_files[1])

    files = sort([f for f in readdir(folder, join=true) if occursin("piv", f)], lt=natural)

    x, y, z, u_dummy = load(files[1], "x", "y", "z", "u")
    x = Int.(x)[1:2:end]
    y = Int.(y)[1:2:end]
    z = Int.(z)[1:2:end]

    u_growth = zeros(Float32, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])
    v_growth = zeros(Float32, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])
    w_growth = zeros(Float32, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])
    u_dispersal = zeros(Float32, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])
    v_dispersal = zeros(Float32, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])
    w_dispersal = zeros(Float32, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])
    flags_tot = zeros(Bool, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])

    for i in 1:8 
        u, v, w, flags = load(files[i], "u", "v", "w", "flags")
        u[flags .> 0] .= 0
        v[flags .> 0] .= 0
        w[flags .> 0] .= 0
        u_growth += permutedims(u, (2,1,3))
        v_growth += permutedims(v, (2,1,3)) 
        w_growth += permutedims(w, (2,1,3))
        flags_tot += permutedims(flags, (2,1,3))
    end
    u_growth[flags_tot .> 3] .= NaN 
    v_growth[flags_tot .> 3] .= NaN  
    w_growth[flags_tot .> 3] .= NaN
    
    for i in 9:length(files)
        u, v, w, flags = load(files[i], "u", "v", "w", "flags")
        u[flags .> 0] .= 0
        v[flags .> 0] .= 0
        w[flags .> 0] .= 0
        u_dispersal += permutedims(u, (2,1,3))
        v_dispersal += permutedims(v, (2,1,3)) 
        w_dispersal += permutedims(w, (2,1,3))
        #flags_tot += permutedims(flags, (2,1,3))
    end
    u_dispersal[flags_tot .> 3] .= NaN 
    v_dispersal[flags_tot .> 3] .= NaN  
    w_dispersal[flags_tot .> 3] .= NaN

    u_growth = u_growth[1:2:end, 1:2:end, 1:2:end] .* 0.065
    v_growth = v_growth[1:2:end, 1:2:end, 1:2:end] .* 0.065
    w_growth = w_growth[1:2:end, 1:2:end, 1:2:end] .* 0.065
    u_dispersal = u_dispersal[1:2:end, 1:2:end, 1:2:end] .* 0.065
    v_dispersal = v_dispersal[1:2:end, 1:2:end, 1:2:end] .* 0.065
    w_dispersal = w_dispersal[1:2:end, 1:2:end, 1:2:end] .* 0.065

    u_growth = mapwindow(median, u_growth, (3,1,1))
    v_growth = mapwindow(median, v_growth, (3,1,1))
    w_growth = mapwindow(median, w_growth, (3,1,1))
    u_dispersal = mapwindow(median, u_dispersal, (3,1,1))
    v_dispersal = mapwindow(median, v_dispersal, (3,1,1))
    w_dispersal = mapwindow(median, w_dispersal, (3,1,1))


    center_of_mass, com_images = compute_com(mask_t1, u_growth)

    radial_component_growth = calculate_radial_component(u_growth, v_growth, w_growth.*4, 8, 8, 8, center_of_mass)
    radial_component_dispersal = calculate_radial_component(u_dispersal, v_dispersal, w_dispersal.*4, 8, 8, 8, center_of_mass)

    fig = Figure(size=(7*72,3*72))
    plot_xlabel = "x (µm)"
    plot_ylabel = "z (µm)"
    ax = CairoMakie.Axis(fig[1,1])
    ax2 = CairoMakie.Axis(fig[1,2])
    ar = arrows!(ax, (x.-1 .- com_images[1]).*0.065, (z.-1).*4 .* 0.065, u_growth[:,15,:], w_growth[:,15,:].*4, arrowsize=4, arrowcolor=vec(radial_component_growth[:,15,:]), linecolor=vec(radial_component_growth[:,15,:]),colorrange=(-2,1), colormap=:berlin)
    ar2 = arrows!(ax2, (x.-1 .- com_images[1]).*0.065, (z.-1).*4 .* 0.065, u_dispersal[:,15,:], w_dispersal[:,15,:].*4, arrowsize=4, arrowcolor=vec(radial_component_dispersal[:,15,:]), linecolor=vec(radial_component_dispersal[:,15,:]),colorrange=(-2,1), colormap=:berlin)
    writedlm("$(plots_folder)/growth_u.csv", u_growth[:,15,:], ",")
    writedlm("$(plots_folder)/growth_v.csv", w_growth[:,15,:], ",")
    writedlm("$(plots_folder)/color_growth.csv", radial_component_growth[:,15,:], ",")
    writedlm("$(plots_folder)/dispersal_u.csv", u_dispersal[:,15,:], ",")
    writedlm("$(plots_folder)/dispersal_v.csv", w_dispersal[:,15,:], ",")
    writedlm("$(plots_folder)/color_dispersal.csv", radial_component_dispersal[:,15,:], ",")
    Colorbar(fig[:, end+1], limits=(-2,1), label="Radial displacement (µm)", colormap=:berlin)
	ax.xlabel = plot_xlabel
    ax.ylabel = plot_ylabel
    ax.title = ""
    ax.rightspinevisible = false
    ax.topspinevisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
	ax2.xlabel = plot_xlabel
    ax2.ylabel = plot_ylabel
    ax2.title = ""
    ax2.rightspinevisible = false
    ax2.topspinevisible = false
    ax2.xgridvisible = false
    ax2.ygridvisible = false
    ax.title = "Growth"
    ax2.title = "Dispersal"
    ylims!(ax, 0, nothing)
    ylims!(ax2, 0, nothing)
    save(plots_folder*"/vector_example.pdf", fig)
end
main()
