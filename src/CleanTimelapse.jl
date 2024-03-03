using Images: channelview, Gray, N0f16, mapwindow 
using NaturalSort: sort, natural
using ImageMorphology: closing, opening, label_components, component_lengths
using StatsBase: median
using FileIO: load, save

include("Preprocessing.jl")
include("RemovePlanktonicCells.jl")
using .Preprocessing
using .RemovePlanktonicCells

function load_raw_images(file_directory, file)
    timeseries = load(file_directory*file)
    height, width, slices, frames = size(timeseries)
    aspect_ratio = 0.3/0.065 
    timeseries = channelview(reshape(timeseries, (height, width, slices, frames)))
    return timeseries, aspect_ratio
end

function time_dependent_processing!(prev_mask, prev_img, curr_mask, curr_img, first_zstack_mask, zstack_mask, t, files, dir, 
                                   intensity_thresholds, height, width, slices, frames, half_window, timepoint_thresh)
    if t < timepoint_thresh || timepoint_thresh < 2
        curr_cv = nothing
        mask_prep_prethresh!(prev_mask, curr_mask, prev_img, curr_img, t,  
                             intensity_thresholds, height, width, slices)
        if t == 2
            @views first_zstack_mask = closing(opening(sum(prev_mask, dims=3)[:,:,1] .> 10, strel_circle(6)), strel_circle(12))
            first_zstack_mask = label_components(first_zstack_mask)
            areas = component_lengths(first_zstack_mask)
            max_label = argmax(areas[1:end])
            first_zstack_mask[first_zstack_mask .!= max_label] .= 0
            first_zstack_mask[first_zstack_mask .== max_label] .= 1 
            prev_mask .*= first_zstack_mask
            prev_img .*= prev_mask
            save("$dir/mask_1.tif", reinterpret(Gray{Bool}, prev_mask))
            save("$dir/noplank_1.tif", prev_img)
        end
    else
        curr_img, curr_mask, curr_cv = mask_prep_postthresh(curr_mask, curr_img, t, files, dir, 
                                                            intensity_thresholds, height, width, 
                                                            slices, frames, half_window)
    end
    isolate_biofilm!(curr_mask, prev_mask, curr_cv, slices)
    if t < timepoint_thresh || timepoint_thresh < 2 
        @views zstack_mask = opening(sum(curr_mask, dims=3)[:,:,1] .> 10, 
                                             strel_circle(6))
        zstack_mask = label_components(zstack_mask)
        areas = component_lengths(zstack_mask)
        max_label = argmax(areas[1:end])
        zstack_mask[zstack_mask .!= max_label] .= 0
        zstack_mask[zstack_mask .== max_label] .= 1 

        zstack_mask = zstack_mask .| first_zstack_mask
        zstack_mask = closing(zstack_mask, strel_circle(12))
        for j in 1:(div(slices, 2))
            curr_mask[:,:,j] .*= zstack_mask
        end
        for j in 1:slices
            curr_mask[:,:,j] = closing(curr_mask[:,:,j], strel_circle(6))
        end
        first_zstack_mask = zstack_mask
    end
end

function processing(timeseries_file, file_directory, dir, cell_threshold, half_window)
    timeseries, aspect_ratio = load_raw_images(file_directory, timeseries_file)
    height, width, slices, frames = size(timeseries)
    
    noback = similar(timeseries)
    rolling_ball!(timeseries, noback, slices, frames)

    timeseries = nothing

    registered = similar(noback)
    center = frames ÷ 2
    register!(noback, registered, frames, center)
    println("Finished registration")
    noback = nothing
    cropped = reinterpret(Gray{N0f16}, crop(registered))
    height, width, slices, frames = size(cropped)
    registered = nothing
    timepoint_thresh = timepoint_threshold(cropped, height, width, slices, frames, cell_threshold)
    
    write_images!(cropped, frames, dir, aspect_ratio)
    intensity_thresholds = Array{Gray{N0f16}, 1}(undef, slices)
    mask_thresholds!(cropped, intensity_thresholds, slices, cell_threshold)
    cropped = nothing

    files = sort([f for f in readdir(dir) if occursin("stack", f)], 
                 lt=natural)

    prev_mask = zeros(Bool, (height, width, slices))
    @views file1 = files[1]
    prev_img = load("$dir/$file1")
    first_zstack_mask = zeros(Bool, (height, width))
    zstack_mask = zeros(Bool, (height, width))
    for t in 2:frames
        @views file = files[t]
        curr_img = load("$dir/$file")
        curr_mask = zeros(Bool, (height, width, slices))
        time_dependent_processing!(prev_mask, prev_img, curr_mask, curr_img, first_zstack_mask, zstack_mask, t, files, dir, 
                                   intensity_thresholds, height, width, slices, frames, half_window, timepoint_thresh)
        curr_img .*= curr_mask 
        prev_mask = curr_mask 
        save("$dir/mask_$(t).tif", reinterpret(Gray{Bool}, curr_mask))
        save("$dir/noplank_$(t).tif", curr_img)
    end
end

function main()
    cell_threshold = 3.0e-5 # Depends on imaging setup!!
    window_size = 5 # For standard deviation in planktonic cell removal
    half_window = window_size ÷ 2

    ##############################
    # Preprocessing
    ##############################

    file_directory = "/mnt/h/Dispersal/test/"
    timeseries_files = [f for f in readdir(file_directory) if occursin("denoised", f)]

    for timeseries_file in timeseries_files
        dir = file_directory*timeseries_file[1:end-16]*"processed"
        if isdir(dir)
            rm(dir; recursive = true)
        end
        mkdir(dir) 
        processing(timeseries_file, file_directory, dir, cell_threshold, half_window)
    end
end
main()
