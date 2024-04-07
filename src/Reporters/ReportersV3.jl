using Images: imfilter, mapwindow, adjust_histogram!, LinearStretching, Gray, N0f16, N0f8,
              Kernel, warp, axes, KernelFactors, RGB 
using HistogramThresholding: find_threshold, Yen, Balanced 
using IntegralArrays: IntegralArray
using IntervalSets: width, leftendpoint, rightendpoint, Interval, ±
using ImageView
using StatsBase: mean, median, quantile
using TiffImages
using FileIO
using NaturalSort: sort, natural
using FLoops

function mean_filter!(X, length_scale)
    iX = IntegralArray(X)
    @inbounds for i in CartesianIndex(1,1):CartesianIndex(size(X))
        x, y = i.I
        x_int = x±length_scale
        y_int = y±length_scale
        x_int = Interval(max(leftendpoint(x_int), 1), 
                         min(rightendpoint(x_int), size(X)[1]))
        y_int = Interval(max(leftendpoint(y_int), 1), 
                         min(rightendpoint(y_int), size(X)[2]))
        X[i] = iX[x_int, y_int]/(width(x_int)*width(y_int))
    end
    return nothing 
end

function normalize_local_contrast(img, blockDiameter, fpMean)
    img_copy = copy(img)
    length_scale = Int((blockDiameter-1)/2)
    mean_filter!(img_copy, length_scale)
    img -= img_copy 
    img .+= fpMean 
    @. img[img < 0.0] = 0.0
    @. img[img > 1.0] = 1.0
    return img 
end

function compute_mask!(stack, masks, ntimepoints, blockDiameter, fpMean, constitutive)
    prev_thresh = 0
    @inbounds for t in 1:ntimepoints
        @views img = imfilter(normalize_local_contrast(stack[:,:,t], blockDiameter, fpMean), Kernel.gaussian(2))
        imshow(img)
        thresh = find_threshold(img, Balanced())
        if thresh < prev_thresh
            thresh = prev_thresh
        end
        mask = img .> thresh
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
    blockDiameter = 401
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
        @views output_stack = images[:,:,constitutive,:]
        masks = zeros(Bool, size(output_stack))
		fpMax = maximum(output_stack) 
		fpMin = minimum(output_stack) 
		fpMean = (fpMax - fpMin) / 2.0 + fpMin
        compute_mask!(output_stack, masks, 
                      nframes, blockDiameter, fpMean, constitutive)
        @views reporter_stack = images[:,:,reporter,:]
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
