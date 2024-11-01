using TiffImages
using LaTeXStrings
using Makie
using GLMakie
using CairoMakie
CairoMakie.activate!(type = "pdf")
set_theme!(fonts=(; bold="TeX Gyre Heros Makie"))
using NaturalSort: sort, natural
using StatsBase
using ColorTypes: Gray, N0f16
using ImageMorphology: label_components, component_lengths, component_centroids
using Images: imresize, distance_transform, feature_transform
using ImageFiltering
using ImageView
using HistogramThresholding: find_threshold, Otsu
using FLoops
using DelimitedFiles
using Interpolations
using FileIO

round_up(x, multiple) = ceil(x / multiple) * multiple

function mask_dispersal_images(downsampled, images_folder)
    mask = zeros(Bool, size(downsampled))
    for i in 1:size(downsampled, 3)
        thresh = find_threshold(downsampled[:,:,i,:], Otsu())
        thresh = max(thresh, 1e-9)
        @views mask[:,:,i,:] = downsampled[:,:,i,:] .> thresh
    end
    for i in 1:size(downsampled, 4)
        TiffImages.save(images_folder*"/downsampled_mask_$(i).tif", Gray.(mask[:,:,:,i]))
    end
    return mask
end

function gaussian_downsample(image_files, first_index, end_index)
    first_img_dummy = TiffImages.load(image_files[first_index])
    first_img_resized = imresize(first_img_dummy, ratio=(1/4, 1/4, 0.3/0.065/4))
    height, width, depth = size(first_img_resized)
    downsampled = Array{Float32, 4}(undef, height, width, depth, end_index-first_index+1)
    for i in first_index:end_index
        img = TiffImages.load(image_files[i])
        img_resized = imresize(img, ratio=(1/4, 1/4, 0.3/0.065/4))
        downsampled[:, :, :, i-first_index+1] = imfilter(img_resized, Kernel.gaussian((3,3,3)))
    end
    return downsampled
end

function radial_averaging(first_index, end_index, bin_interval, dispersal_mask)
    # Step 1: get the center at the substrate of the biofilm at the first timepoint
    first_mask = dispersal_mask[:,:,:,1] 
    labels = label_components(first_mask)
    volumes = component_lengths(labels)
    centers = component_centroids(labels)
    # Biofilm = largest connected component
    center = centers[argmax(volumes[1:end])]
    center = (round(Int, center[1]), round(Int, center[2]), 1)
    # Step 2: get the distance of each voxel from the center
    center_distance = zeros(Bool, size(labels))
    center_distance[center[1], center[2], center[3]] = true
    center_distance = distance_transform(feature_transform(center_distance))
    # Get the maximum distance from the center
    max_distance = round_up(maximum(center_distance), bin_interval)
    bins = 0:bin_interval:max_distance
    nbins = length(bins)
    ntimepoints = end_index - first_index
    # Initialize kymograph (transposed) and boundary array
    data_matrix = zeros(nbins-1, ntimepoints)
    boundary = zeros(ntimepoints)
    # Step 3: for each timepoint, the density in a given bin is the fraction of 1's inside the bin
    # boundary is the maximum distance from the center at which a 1 is present in the mask (units of voxels)
    for t in 1:ntimepoints
        @views curr_mask = dispersal_mask[:,:,:,t]
        for i in 1:size(data_matrix, 1) 
            # Find the voxels inside the bin by using the distance transform from the center
            # Calculate the fraction of 1's in the bin by indexing the mask with the voxels inside the bin
            data_matrix[i, t] = mean(curr_mask[findall(x -> bins[i] <= x <= bins[i+1], center_distance)])
        end
        boundary[t] = maximum(center_distance .* curr_mask)
    end
    return data_matrix, boundary
end

function main()
    plot_xlabel = "Time (h)"
    plot_ylabel = "Distance from center (µm)"
    plot_title = "Data"

    master_directory = "/mnt/h/Dispersal"
    image_folders = filter(isdir, readdir(master_directory, join=true))
    image_folders = [f for f in image_folders if !occursin("rbmA", f)]
    filter!(folder->folder≠master_directory*"/Plots", image_folders)
    plots_folder = "/mnt/h/Dispersal/Plots"

    for images_folder in image_folders
        plot_filename = basename(images_folder)*"_data_intensity" 
        displacements_folder = "$(images_folder)/Displacements"
        image_files = sort([f for f in readdir(images_folder, join=true) if occursin("noplank", f)], 
                                 lt=natural)
        ntimepoints = length(image_files)
        net = readdlm(plots_folder*"/"*basename(images_folder)*".csv", ',', Int)[1:end,1]
        # Get first and last indices for dispersal
        first_index = argmax(net)
        end_index = min(first_index + 45, ntimepoints)
        # Downsample the relevant images
        downsampled = gaussian_downsample(image_files, first_index, end_index)
        # Mask the downsampled images
        dispersal_mask = mask_dispersal_images(downsampled, images_folder)
        # Radially average the result to generate the kymograph (units of density)
        # Chosen bin side length = 8 pixels -> 0.065*4 µm / pixel -> 2.08 µm bin side length
        data_matrix, boundary = radial_averaging(first_index, end_index, 8, dispersal_mask)
        # Subtract densities at consecutive timepoints to get density change / 10 min (*6 to get /hr) 
        data_matrix = diff(data_matrix, dims=2) .* 6
        boundary = boundary ./ 8 # Convert boundary from units of voxels to units of bins
        data_matrix[data_matrix .> 0] .= 0 # Isolate dispersing regions
        writedlm("$(plots_folder)/$(plot_filename).csv", data_matrix, ",")
        n = 10 # Desired ytick interval in µm
        ytick_interval = n/(0.065*4*8) # ytick interval =  desired µm interval / bin side length in µm
        xs = 0:6:size(data_matrix, 2)-1 # xticks every hour
        ys = 0:ytick_interval:size(data_matrix, 1)-1
        fig = Figure(size=(5*72, 3*72))
        ax = Axis(fig[1, 1])
        colormap = cgrad(["#cf34eb", :white])
        hm = heatmap!(ax, 0:size(data_matrix,2), 0:size(data_matrix,1), 
                      transpose(data_matrix), colormap=colormap, padding=(0.0, 0.0))
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
