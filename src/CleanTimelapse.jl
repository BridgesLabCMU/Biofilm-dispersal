using Images: channelview, Gray, N0f16, mapwindow 
using NaturalSort: sort, natural
using ImageMorphology: label_components, component_lengths, closing
using StatsBase: median
using FileIO: load, save
using HistogramThresholding: find_threshold, Otsu
using DelimitedFiles
using PythonCall
morph = pyimport("skimage.morphology")

include("Preprocessing.jl")
include("RemovePlanktonicCells.jl")
using .Preprocessing
using .RemovePlanktonicCells

function load_raw_images(file_directory, file)
    timeseries = load(file_directory*file)
    height, width, slices, frames = size(timeseries)
    timeseries = channelview(reshape(timeseries, (height, width, slices, frames)))
    return timeseries
end

function time_dependent_processing(prev_mask, prev_img, curr_mask, curr_img, first_zstack_mask, zstack_mask, t, files, dir, 
                                   intensity_thresholds, height, width, slices, frames, half_window, timepoint_thresh, zstack_thresh)
    closed_prev = nothing
    if t < timepoint_thresh || timepoint_thresh < 2
        mask_prep!(prev_mask, curr_mask, prev_img, curr_img, t,  
                             intensity_thresholds, height, width, slices)
        if t == 2
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
        closed_prev = pyconvert(Array, morph.isotropic_closing(prev_mask, radius=25))
        mask_prep!(prev_mask, curr_mask, prev_img, curr_img, t,  
                             intensity_thresholds, height, width, slices)
    end
    isolate_biofilm!(curr_mask, prev_mask, closed_prev, slices, t, timepoint_thresh)
    if t < timepoint_thresh || timepoint_thresh < 2 
        @views zstack = sum(curr_mask, dims=3)[:,:,1]
        @views zstack_mask = pyconvert(Array, morph.isotropic_opening(zstack .> zstack_thresh/2, radius=6))
        zstack_mask = label_components(zstack_mask)
        areas = component_lengths(zstack_mask)
        max_label = argmax(areas[1:end])
        zstack_mask[zstack_mask .!= max_label] .= 0
        zstack_mask[zstack_mask .== max_label] .= 1 

        zstack_mask = zstack_mask .| first_zstack_mask
        zstack_mask = pyconvert(Array, morph.isotropic_closing(zstack_mask, radius=12))
        for j in 1:(div(slices, 2))
            curr_mask[:,:,j] .*= zstack_mask
        end
        for j in 1:slices
            curr_mask[:,:,j] = pyconvert(Array, morph.isotropic_closing(curr_mask[:,:,j], radius=6))
        end
    else
        for j in 1:(div(slices, 2))
            curr_mask[:,:,j] .*= zstack_mask
        end
        for j in 1:slices
            curr_mask[:,:,j] = pyconvert(Array, morph.isotropic_closing(curr_mask[:,:,j], radius=6))
        end
    end
    return zstack_mask, first_zstack_mask, curr_mask, curr_img
end

function processing(timeseries_file, file_directory, dir, cell_threshold, half_window, registered_flag, frames)

    if !registered_flag
        println("This didn't work")
        timeseries = load_raw_images(file_directory, timeseries_file)
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
        writedlm("$(dir)/timepoint_thresh.csv", [timepoint_thresh], ",")
        
        write_images!(cropped, frames, dir)
        intensity_thresholds = Array{Gray{N0f16}, 1}(undef, slices)
        mask_thresholds!(cropped, intensity_thresholds, slices, cell_threshold)
        intensity_thresholds = Float64.(intensity_thresholds)
        writedlm("$(dir)/intensity_thresholds.csv", intensity_thresholds, ",")
        cropped = nothing
    end

    files = sort([f for f in readdir(dir) if occursin("stack", f)], 
                 lt=natural)
    @views intensity_thresholds = readdlm("$(dir)/intensity_thresholds.csv", ',', Float64)[1:end,1]
    @views timepoint_thresh = readdlm("$(dir)/timepoint_thresh.csv", ',', Float64)[1,1]
    @views file1 = files[1]
    prev_img = load("$dir/$file1")
    height, width, slices = size(prev_img)
    prev_mask = zeros(Bool, (height, width, slices))
    for i in 1:slices
        prev_mask[:,:,i] = closing(prev_img[:,:,i] .> intensity_thresholds[i], 
                                          strel_circle(6))
    end
    @views first_zstack = sum(prev_mask, dims=3)[:,:,1]
    zstack_thresh = 15 #find_threshold(first_zstack, Otsu())
    @views first_zstack_mask = pyconvert(Array, morph.isotropic_closing(pyconvert(Array, 
                                                morph.isotropic_opening(first_zstack .> zstack_thresh/2, 
                                                        radius=6)), radius=12))
    zstack_mask = zeros(Bool, (height, width))
    for t in 2:frames
        @views file = files[t]
        curr_img = load("$dir/$file")
        curr_mask = zeros(Bool, (height, width, slices))
        zstack_mask, first_zstack_mask, curr_mask, curr_img = time_dependent_processing(prev_mask, prev_img, 
                                                                                         curr_mask, curr_img, 
                                                                                         first_zstack_mask, zstack_mask, t, files, 
                                                                                         dir, intensity_thresholds, height, width,
                                                                                         slices, frames, half_window,
                                                                                         timepoint_thresh, zstack_thresh)
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
        if !isdir(dir)
            mkdir(dir) 
        end
        registered_files = [f for f in readdir(dir) if occursin("stack", f)]
        if isempty(registered_files)
            registered_flag = false
            frames = nothing
        else
            registered_flag = true
            frames = length(registered_files)
        end
        processing(timeseries_file, file_directory, dir, cell_threshold, half_window, registered_flag, frames)
    end
end
main()
