using Makie
using GLMakie
using CairoMakie
CairoMakie.activate!(type = "pdf")
set_theme!(fonts=(; bold="TeX Gyre Heros Makie"))
using LaTeXStrings
using TiffImages
using NaturalSort: sort, natural
using StatsBase
using ColorTypes: Gray, N0f16
using ImageMorphology: label_components, component_lengths, component_centroids
using Images: imresize, distance_transform, feature_transform
using HistogramThresholding: find_threshold, Otsu
using FLoops
using DelimitedFiles

round_up(x, multiple) = ceil(x / multiple) * multiple

function load_images(image_files, first_index, end_index)
    dummy_image = TiffImages.load(image_files[first_index]; lazyio=true)
    height, width, depth = size(dummy_image)
    images = Array{Bool, 4}(undef, height, width, depth, end_index-first_index+1)
    for t in first_index:end_index
        images[:,:,:,t-first_index+1] = TiffImages.load(image_files[t-first_index+1]) .> 0
    end
    return images
end

function radial_averaging(first_index, end_index, bin_interval, masks)
    first_image = masks[:,:,:,1] 
    labels = label_components(first_image)
    volumes = component_lengths(labels)
    centers = component_centroids(labels)
    center = centers[argmax(volumes[1:end])]
    center = (round(Int, center[1]), round(Int, center[2]), 1)
    center_distance = zeros(Bool, size(labels))
    center_distance[center[1], center[2], center[3]] = true
    center_distance = distance_transform(feature_transform(center_distance))
    max_distance = round_up(maximum(center_distance), bin_interval)
    bins = 0:bin_interval:max_distance
    nbins = length(bins)
    ntimepoints = end_index - first_index
    data_matrix = zeros(nbins-1, ntimepoints-1)
    boundary = zeros(ntimepoints)
    budget = zeros(Int, ntimepoints-1)
    dispersal_config = zeros(Bool, size(center_distance))
    for t in 1:ntimepoints-1
        @views curr_mask = masks[:,:,:,t]
        @views next_mask = masks[:,:,:,t+1]
        for i in 1:size(data_matrix, 1) 
            curr_density = curr_mask[findall(x -> bins[i] <= x <= bins[i+1], center_distance)]
            next_density = next_mask[findall(x -> bins[i] <= x <= bins[i+1], center_distance)]
            if sum(curr_density) > sum(next_density)
                budget[t] -= (sum(next_density) - sum(curr_density))
            end
        end
    end
    for t in 1:ntimepoints-1
        biofilm_configuration = masks[:,:,:,t]
        dispersal_config .= 0 
        if budget[t] <= sum(biofilm_configuration)
            N_voxels = sample(findall(!iszero, biofilm_configuration), budget[t], replace=false)
            boundary[t] = maximum(center_distance .* (biofilm_configuration .> 0))
            dispersal_config[N_voxels] .= 1 
            for i in 1:size(data_matrix, 1) 
                data_matrix[i, t] = -1*mean(dispersal_config[findall(x -> bins[i] <= x <= bins[i+1], center_distance)])
            end
        end
    end
    found_zero = false
    for i in eachindex(boundary)
        if found_zero
            boundary[i] = 0 
        elseif boundary[i] == 0 
            found_zero = true
        end
	end
    return data_matrix, boundary
end

function main()
    plot_xlabel = "Time (h)"
    plot_ylabel = "Distance from center (µm)"
    plot_title = "Random model"

    master_directory = "/mnt/h/Dispersal"
    image_folders = filter(isdir, readdir(master_directory, join=true))
    image_folders = [f for f in image_folders if !occursin("Plots", f)]
    image_folders = [f for f in image_folders if occursin("rbmA", f)]
    filter!(folder->folder≠master_directory*"/Plots", image_folders)
    plots_folder = "/mnt/h/Dispersal/Plots"

    for images_folder in image_folders
        plot_filename = basename(images_folder)*"_random_net_downsampled" 
        files = sort([f for f in readdir(images_folder, join=true) if occursin("downsampled_mask", f)], 
                                 lt=natural)
        all_files = sort([f for f in readdir(images_folder, join=true) if occursin("stack", f)], 
                                 lt=natural)
        ntimepoints = length(all_files)
        if occursin("rbmA", images_folder)
            net = readdlm(plots_folder*"/"*basename(images_folder)*".csv", ',', Int)[1:end,1] .+ 0.01
            data_max = maximum(net)
            fc = net ./ data_max
            min_fc = minimum(fc)
            fc = fc .- min_fc
            end_index = findfirst(x->x<0.001, fc) 
            first_index = end_index - 45 
        else
            net = readdlm(plots_folder*"/"*basename(images_folder)*".csv", ',', Int)[1:end,1]
            first_index = argmax(net)
            end_index = min(first_index + 45, ntimepoints)
        end
        masks = load_images(files, first_index, end_index)
        data_matrix, boundary = radial_averaging(first_index, end_index, 8, masks)
        data_matrix .*= 6
        boundary ./= 8
        writedlm("$(plots_folder)/$(plot_filename).csv", data_matrix, ",")
        n = 10 
        ytick_interval = n/(0.065*4*8) # interval in µm / bin side length in µm
        xs = 0:6:size(data_matrix, 2)-1
        ys = 0:ytick_interval:size(data_matrix, 1)-1
        fig = Figure(size=(5*72, 3*72))
        ax = Axis(fig[1, 1])
        colormap = cgrad(["#cf34eb", :white])
        vmin = quantile(vec(data_matrix), 0.03)
        hm = heatmap!(ax, 0:size(data_matrix,2), 0:size(data_matrix,1), 
                      transpose(data_matrix), colormap=colormap, colorrange=(vmin, 0), padding=(0.0, 0.0))
        lines!(ax, 0:size(data_matrix, 2), boundary, color=:black)
        Colorbar(fig[:, end+1], hm, label="Density change/h")
        ax.xticks = xs
        ax.yticks = ys
        ax.xtickformat=values->string.([Int(div(v,6)) for v in values])
        ax.ytickformat=values->string.([Int(v*n/ytick_interval) for v in values])
        ax.xlabel = plot_xlabel
        ax.ylabel = plot_ylabel
        ax.title = plot_title
        ax.xgridvisible = false
        ax.ygridvisible = false
        ylims!(ax, 0, maximum(boundary)+2)
        xlims!(ax, 0, size(data_matrix, 2)-1)
        save("$plots_folder/$plot_filename"*".pdf", fig)
    end
end

main()

