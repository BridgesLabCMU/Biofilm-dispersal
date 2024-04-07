using Images: imfilter, mapwindow, adjust_histogram!, LinearStretching, Gray, N0f16, N0f8,
              Kernel, warp, axes, KernelFactors, RGB 
using HistogramThresholding: find_threshold, Yen 
using ImageView
using StatsBase: mean, median, quantile
using TiffImages
using FileIO
using NaturalSort: sort, natural
using FLoops
using PythonCall

morph = pyimport("skimage.morphology")

function compute_mask!(stack, masks, ntimepoints, constitutive)
    @inbounds for t in 1:ntimepoints
        @views img = imfilter(stack[:,:,t], Kernel.gaussian(2)) 
        mask = img .> find_threshold(img, Yen()) 
        #mask = pyconvert(Array, morph.isotropic_opening(mask, radius=5))
		masks[:,:,t] = mask
	end
end

function output_images!(stack, masks, overlay, output_dir, replicate)
    imshow(masks)
    stack .*= masks
    flat_stack = vec(stack)
    img_min = quantile(flat_stack, 0.0035)
    img_max = quantile(flat_stack, 0.9965)
    adjust_histogram!(stack, LinearStretching(src_minval=img_min, src_maxval=img_max, 
                                              dst_minval=0, dst_maxval=1))
    TiffImages.save("$output_dir/constitutive_"*replicate*".tif", stack)
    @inbounds for i in CartesianIndices(stack)
        gray_val = RGB{N0f8}(stack[i], stack[i], stack[i])
        overlay[i] = masks[i] ? RGB{N0f8}(1, 0, 0) : gray_val
    end
    TiffImages.save("$output_dir/constitutive_mask_"*replicate*".tif", overlay)
end

function main()
    dir = "/mnt/f/Sandhya_Imaging/Time_Lapses/Reporters/vpsL/" # Folder containing the images for strain X 
    files = sort([f for f in readdir(dir, join=true) if occursin(".ome.tif", f)], lt=natural)
    shift_thresh = 100 
    sig = 2 
    constitutive = 2 # index for constitutive channel
    reporter = 1 # index for reporter channel
    ntimepoints = 17
    if isdir("$dir/results_images")
        rm("$dir/results_images"; recursive = true)
    end
    if isdir("$dir/results_data")
        rm("$dir/results_data"; recursive = true)
    end
    output_img_dir = "$dir/results_images"
    output_data_dir = "$dir/results_data"
    mkdir(output_img_dir)
    mkdir(output_data_dir)
    output_file = "$output_data_dir/data.jld2"
    data_matrix = Array{Float64, 2}(undef, ntimepoints, length(files))
    for (j, file) in enumerate(files)
        height, width, channels, nslices, nframes = size(FileIO.load(file))
        if nframes != ntimepoints
            error("Number of timepoints in the image stack ($nframes) does not match the expected number of timepoints ($ntimepoints)")
        end
        @views images = Float64.(sum(reshape(TiffImages.load(file), 
                                             (height, width, nslices, channels, nframes)), dims=3)[:,:,1,:,:])
        output_stack = images[:,:,constitutive,:]
        @views first_slice = output_stack[:,:,1]
        output_stack = output_stack .- first_slice[:,:,:]
        output_stack[output_stack .< 0] .= 0
        @views masks = zeros(Bool, size(output_stack))
        compute_mask!(output_stack, masks, 
                      nframes, constitutive)
        @views first_slice = images[:,:,reporter,1]
        reporter_stack = images[:,:,reporter,:] .- first_slice[:,:,:]
        reporter_stack[reporter_stack .< 0] .= 0
        let output_stack = output_stack 
            @floop for t in 1:nframes
                @inbounds signal = @views mean((reporter_stack[:,:,t] ./ output_stack[:,:,t]) .* masks[:,:,t])
                @inbounds data_matrix[t, j] = signal 
            end
        end
        output_stack = Gray{N0f8}.(output_stack)
        overlay = zeros(RGB{N0f8}, size(output_stack)...)
        output_images!(output_stack, masks, overlay, output_img_dir, string(j))
    end
    data_matrix .= ifelse.(isnan.(data_matrix), 0, data_matrix)
    FileIO.save(output_file, Dict("data" => data_matrix))
end
main()
