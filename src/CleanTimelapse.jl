using TiffImages: load, save, ifds, IMAGEDESCRIPTION
using Images: channelview, Gray, N0f16 
using NaturalSort: sort, natural
using StatsBase: mean, std
using ImageMorphology: closing, tophat!

include("Preprocessing.jl")
include("RemovePlanktonicCells.jl")
using .Preprocessing
using .RemovePlanktonicCells

function load_raw_images(file_directory, file)
    timeseries = load(file_directory*file)
    height, width, zt = size(timeseries)
    ImageJ_metadata = first(ifds(timeseries))[IMAGEDESCRIPTION].data
    xy_res = 0.065
    z_spacing = findfirst("spacing=", ImageJ_metadata)
    nloop = findfirst("loop=", ImageJ_metadata)
    z_res = tryparse(Float32, ImageJ_metadata[(z_spacing[8]+1):(nloop[1]-1)])
    aspect_ratio = z_res / xy_res
    slice_loc = findfirst("slices=", ImageJ_metadata)
    frames_loc = findfirst("frames=", ImageJ_metadata)
    hyperstack_loc = findfirst("hyperstack=", ImageJ_metadata)
    slices = tryparse(Int, ImageJ_metadata[(slice_loc[7]+1):(frames_loc[1]-1)])
    frames = tryparse(Int, ImageJ_metadata[(frames_loc[7]+1):(hyperstack_loc[1]-1)])
    timeseries = channelview(reshape(timeseries, (height, width, slices, frames)))
    return timeseries, aspect_ratio
end

function main()
    cell_threshold = 3.0e-5 # Depends on imaging setup!!
    window_size = 5 # For standard deviation in planktonic cell removal
    half_window = window_size ÷ 2

    ##############################
    # Preprocessing
    ##############################

    timeseries_file = "cheY_replicate1_denoised.tif" 
    file_directory = "/mnt/h/Dispersal/"
	if isdir(file_directory*timeseries_file[1:end-12]*"processed")
        rm(file_directory*timeseries_file[1:end-12]*"processed"; recursive = true)
	end
    dir = file_directory*timeseries_file[1:end-12]*"processed"
	mkdir(dir) 
    timeseries, aspect_ratio = load_raw_images(file_directory, timeseries_file)
    height, width, slices, frames = size(timeseries)
    registered = similar(timeseries)
    center = frames ÷ 2
    register!(timeseries, registered, frames, center)
    timeseries = nothing
    cropped = reinterpret(Gray{N0f16}, crop(registered))
    height, width, slices, frames = size(cropped)
    @show height, width, slices, frames
    registered = nothing
    timepoint_thresh = timepoint_threshold(cropped, height, width, slices, frames, cell_threshold)
    
    noback = similar(cropped)
    for t in 1:frames
        for i in 1:slices
            @views buffer = similar(cropped[:,:,i,t])
            @views tophat!(noback[:,:,i,t], cropped[:,:,i,t], buffer; r=15)
            buffer = nothing
        end
    end

    cropped = nothing
    write_images!(noback, frames, dir, aspect_ratio)
    intensity_thresholds = Array{Gray{N0f16}, 1}(undef, slices)
    mask_thresholds!(noback, intensity_thresholds, slices, cell_threshold)
    noback = nothing

    ##############################
    # Planktonic cell removal
    ##############################

    files = sort([f for f in readdir(dir) if occursin("stack", f)], 
                 lt=natural)

    prev_mask = zeros(Bool, (height, width, slices))
    for t in 2:frames
        # Preparing to run isolate_biofilm function
        if t < timepoint_thresh || timepoint_thresh < 2
            curr_cv = nothing
            if t == 2
                @views file1 = files[1]
                @views file2 = files[2]
                prev_img = load("$dir/$file1")
                curr_img = load("$dir/$file2")
                curr_mask = zeros(Bool, (height, width, slices))
                for i in 1:slices
                    prev_mask[:,:,i] = closing(prev_img[:,:,i] .> intensity_thresholds[i], 
                                                      strel_circle(6))
                    curr_mask[:,:,i] = curr_img[:,:,i] .> intensity_thresholds[i]
                end
                prev_img .*= prev_mask
                save("$dir/mask_1.tif", reinterpret(Gray{Bool}, prev_mask))
                save("$dir/noplank_1.tif", prev_img)
            else
                @views file = files[t]
                curr_img = load("$dir/$file")
                curr_mask = zeros(Bool, (height, width, slices))
                for i in 1:slices
                    curr_mask[:,:,i] = curr_img[:,:,i] .> intensity_thresholds[i]
                end
            end
        else
            first_index = max(1, t-half_window)
            last_index = min(frames, t+half_window)
            window = Array{N0f16, 4}(undef, height, width, slices, last_index-first_index+1) 
            window_mask = zeros(Bool, (height, width, slices, last_index-first_index+1))
            @views file = files[t]
            curr_img = load("$dir/$file")
            curr_mask = zeros(Bool, (height, width, slices))
           
            for i in 1:slices
                curr_mask[:,:,i] = curr_img[:,:,i] .> intensity_thresholds[i]
            end
            
            for τ in first_index:last_index
                @views file = files[τ]
                window[:,:,:,τ-first_index+1] = load("$dir/$file")
            end
            
            for τ in size(window, 4)
                for i in 1:slices
                    window_mask[:,:,i,τ] = window[:,:,i,τ] .> intensity_thresholds[i]
                end
            end

            window = nothing
            mask_float = Float32.(window_mask)
            window_mask = nothing
            mask_mean = mean(mask_float, dims=4) 
            mask_std = std(mask_float, dims=4) 
            @views curr_cv = @. ifelse(iszero(mask_mean), mask_std/mask_mean, 0.0)[:,:,:,1]
            mask_float = nothing
            mask_mean = nothing
            mask_std = nothing
        end

        isolate_biofilm!(curr_img, curr_mask, prev_mask, curr_cv, slices)
        prev_mask = curr_mask 
        save("$dir/mask_$(t).tif", reinterpret(Gray{Bool}, curr_mask))
        save("$dir/noplank_$(t).tif", curr_img)
    end
end
main()
