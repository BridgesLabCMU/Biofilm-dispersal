using TiffImages: load, save
using NaturalSort: sort, natural
using StatsBase
using ColorTypes: Gray, N0f16
using ImageMorphology: label_components, component_lengths, component_centroids
using Images: imresize, WarpedView, distance_transform, feature_transform
using HistogramThresholding: find_threshold, Otsu
using FLoops
using DelimitedFiles
using LinearAlgebra
using CoordinateTransformations

function read_images!(directory, ntimepoints, arr, files)
    @inbounds for t in 1:ntimepoints
        @views file = files[t]
        arr[:, :, :, t] = load("$directory/$file") .> 0
    end
    return nothing 
end

function write_images!(binary_timeseries, frames, dir)
    @floop for t in 1:frames
        @views isotropic = reinterpret(Gray{Bool}, binary_timeseries[:,:,:,t])
        save("$dir/mask_isotropic_$(t).tif", isotropic)
    end
end

function mask_thresholds!(timeseries, intensity_thresholds, slices, cell_threshold)
    @floop for i in 1:slices 
        @views slice_ = timeseries[:, :, i, :]
        otsu_thresh = find_threshold(slice_, Otsu())
        if otsu_thresh/3 < cell_threshold 
            intensity_thresholds[i] = cell_threshold 
        else
            intensity_thresholds[i] = otsu_thresh/3
        end
    end
    return nothing 
end
function main()
    master_directory = "/mnt/h/Dispersal"
    image_folders = filter(isdir, readdir(master_directory, join=true))
    filter!(folder->folder≠master_directory*"/Plots", image_folders)
    image_folders = [f for f in image_folders if occursin("rbmA", f) && occursin("replicate3", f)]
    plots_folder = "/mnt/h/Dispersal/Plots"
    for images_folder in image_folders
        filename = basename(images_folder)
        files = sort([f for f in readdir(images_folder) if occursin("noplank", f)], 
                                 lt=natural)
        ntimepoints = length(files)
        dummy_image = load("$images_folder/$(files[1])"; lazyio=true)
        height, width, slices = size(dummy_image)
        timeseries = zeros(Gray{N0f16}, height, width, slices, ntimepoints)
        read_images!(images_folder, ntimepoints, timeseries, files)
        # Make the timeseries isotropic in z
        timeseries = imresize(timeseries, ratio=(1,1,0.3/0.065,1))
        height, width, slices, ntimepoints = size(timeseries)
        # Separate intensity threshold at each z-slice
        intensity_thresholds = zeros(slices)
        # Find the thresholds
        mask_thresholds!(timeseries, intensity_thresholds, slices, 3.0e-5)
        masks = zeros(Bool, height, width, slices, ntimepoints)
        for i in 1:slices
            masks[:,:,i,:] = @views timeseries[:,:,i,:] .> intensity_thresholds[i]
        end
        timeseries = nothing
        # Write the isotropic masks
        write_images!(masks, ntimepoints, images_folder)
        # Calculate biovolume (sum of all binary voxels for each timepoint) 
        @views net = sum(masks, dims=(1:3))[1,1,1,:] 
        # Write intensity thresholds and biovolumes
        writedlm("$(plots_folder)/$(filename).csv", net, ",")
        writedlm("$(images_folder)/isotropic_intensity_thresholds.csv", intensity_thresholds, ",")
        masks = nothing
    end
end

main()
