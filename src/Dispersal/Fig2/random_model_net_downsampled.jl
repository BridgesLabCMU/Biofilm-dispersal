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

function N_smallest(M, center_distance, biofilm_configuration, n)
    M[center_distance .== 0] .= biofilm_configuration[center_distance .== 0]  
    v = vec(M)
    nonzero_indices = findall(x -> x != 0, v)
    n = min(n, length(nonzero_indices))
    if n == 0
        return CartesianIndex[]
    end
    perm = partialsortperm(v[nonzero_indices], 1:n)
    indices = CartesianIndices(M)[nonzero_indices[perm]]
    return indices
end

function radial_averaging(images_folder, files, first_index, end_index, bin_interval, budget)
    first_image = TiffImages.load(files[first_index])
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
    data_matrix = zeros(nbins-1, ntimepoints)
    biofilm_configuration = TiffImages.load(files[1])
    dispersal_config = zeros(Bool, size(biofilm_configuration))
    boundary = zeros(ntimepoints)
    for t in 1:ntimepoints
        dispersal_config .= 0 
        N_voxels = sample(findall(!iszero, biofilm_configuration), budget[t], replace=false) 
        boundary[t] = maximum(center_distance .* (biofilm_configuration .> 0))
        biofilm_configuration[N_voxels] .= 0 
        dispersal_config[N_voxels] .= 1 
        for i in 1:size(data_matrix, 1) 
            data_matrix[i, t] = -1*mean(dispersal_config[findall(x -> bins[i] <= x <= bins[i+1], center_distance)])
        end
        if t == (ntimepoints-1)
            TiffImages.save(images_folder*"/random_mask_final.tif", Gray{Bool}.(biofilm_configuration.>0))
        end
    end
    return data_matrix, boundary
end

function main()
    plot_xlabel = "Time (h)"
    plot_ylabel = "Distance from center \n (µm)"
    plot_title = "Random model"

    master_directory = "/mnt/h/Dispersal"
    image_folders = filter(isdir, readdir(master_directory, join=true))
    image_folders = [f for f in image_folders if occursin("WT", f)]
    filter!(folder->folder≠master_directory*"/Plots", image_folders)
    plots_folder = "/mnt/h/Dispersal/Plots"

    for images_folder in image_folders
        plot_filename = basename(images_folder)*"_random_net_downsampled" 
        if isfile("$plots_folder/$plot_filename"*".pdf")
            continue
        end
        files = sort([f for f in readdir(images_folder, join=true) if occursin("downsampled_mask", f)], 
                                 lt=natural)
        all_files = sort([f for f in readdir(images_folder, join=true) if occursin("stack", f)], 
                                 lt=natural)
        ntimepoints = length(all_files)
        net = readdlm(plots_folder*"/"*basename(images_folder)*".csv", ',', Int)[1:end,1]
        first_index = argmax(net)
        end_index = min(first_index + 45, ntimepoints)
        masks = load_images(files, first_index, end_index)
        @views budget = diff(sum(masks, dims=(1,2,3))[1,1,1,:]) .* -1
        pos_avg = sum(budget[budget .< 0])[1]
        if !isnan(pos_avg)
            budget[budget .> 0] .+= (round(Int,pos_avg/length(budget[budget .> 0]))) 
        end
        budget[budget .< 0] .= 0
        data_matrix, boundary = radial_averaging(images_folder, files, first_index, end_index, 8, budget)
        data_matrix .*= 6
        data_matrix = data_matrix[:, 1:end-1]
        boundary ./= 8
        writedlm("$(plots_folder)/$(plot_filename).csv", data_matrix, ",")
        n = 10 
        ytick_interval = n/0.065/30
        xs = 0:6:size(data_matrix, 2)-1
        ys = 0:ytick_interval:size(data_matrix, 1)-1
        fig = Figure(size=(5*72, 3*72))
        ax = Axis(fig[1, 1])
        hm = heatmap!(ax, 0:size(data_matrix,2), 0:size(data_matrix,1), 
                      transpose(data_matrix), colormap=:devon, padding=(0.0, 0.0))
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
