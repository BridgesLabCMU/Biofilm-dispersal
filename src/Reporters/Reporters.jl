using Images: imfilter, mapwindow, adjust_histogram!, LinearStretching, Gray, N0f16, N0f8,
              Kernel, warp, axes, KernelFactors, RGB 
using ImageMorphology: label_components, component_lengths
using StatsBase: mean, median, quantile
using FileIO
using NaturalSort: sort, natural
using DataFrames: DataFrame
using CSV: write
using JSON: parsefile
using IntegralArrays: IntegralArray
using CoordinateTransformations: Translation
using IntervalSets: width, leftendpoint, rightendpoint, Interval, ±
using SubpixelRegistration: phase_offset
using FLoops

function thresh(μ, t₀)
    if μ < t₀/1.5 
        return (t₀+0.03) + ((t₀-0.03)-(t₀+0.03))*μ
    else
        return (t₀+0.06) + ((t₀-0.06)-(t₀+0.06))*μ
    end
end

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

function normalize_local_contrast(img, img_copy, blockDiameter, fpMean)
    length_scale = Int((blockDiameter-1)/2)
    mean_filter!(img_copy, length_scale)
    img -= img_copy 
    img .+= fpMean 
    @. img[img < 0.0] = 0.0
    @. img[img > 1.0] = 1.0
    return img 
end

function crop(img_stack)
    @views mask = .!any(isnan.(img_stack), dims=3)[:,:,1]
    @views mask_i = any(mask, dims=2)[:,1]
    @views mask_j = any(mask, dims=1)[1,:]
    i1 = findfirst(mask_i)
    i2 = findlast(mask_i)
    j1 = findfirst(mask_j)
    j2 = findlast(mask_j)
    cropped_stack = img_stack[i1:i2, j1:j2, :]
    return cropped_stack, (i1, i2, j1, j2)
end

function stack_preprocess(img_stack, normalized_stack, registered_stack,
                        blockDiameter, shift_thresh, fpMean, nframes, mxshift, sig)       
    shifts = (0.0, 0.0) 
    @inbounds for t in 1:nframes
        img = img_stack[:,:,t]
        img_copy = img_stack[:,:,t] 
        img_normalized = normalize_local_contrast(img, img_copy, 
                                    blockDiameter, fpMean)
        normalized_stack[:,:,t] = imfilter(img_normalized, Kernel.gaussian(sig))
        if t == 1
            registered_stack[:,:,t] = normalized_stack[:,:,t]
        else
            moving = normalized_stack[:,:,t]
            fixed = normalized_stack[:,:,t-1]
            shift, _, _ = phase_offset(fixed, moving, upsample_factor=1)
            if sqrt(shift[1]^2 + shift[2]^2) >= mxshift
                shift = Translation(shifts[1], shifts[2])
                registered_stack[:,:,t] = warp(moving, shift, axes(fixed))
                img_stack[:,:,t] = warp(img_stack[:,:,t], shift, axes(fixed))
            else
                shift = Tuple([-1*shift[1], -1*shift[2]])
                shift = shift .+ shifts
                shifts = shift
                shift = Translation(shift[1], shift[2])
                registered_stack[:,:,t] = warp(moving, shift, axes(fixed))
                img_stack[:,:,t] = warp(img_stack[:,:,t], shift, axes(fixed))
            end
        end
    end
    processed_stack, crop_indices = crop(registered_stack)
    row_min, row_max, col_min, col_max = crop_indices
    img_stack = img_stack[row_min:row_max, col_min:col_max, :]
    return img_stack, processed_stack
end

function compute_mask!(stack, raw_stack, masks, 
                        sig, fixed_thresh, ntimepoints)
    @inbounds for t in 1:ntimepoints
        @views img = stack[:,:,t]
        plank_mask = img .> fixed_thresh 
        @views plank = plank_mask .* raw_stack[:,:,t] 
        flattened_plank = vec(plank)
        plank_pixels = filter(x -> x != 0, flattened_plank)
        plank_avg = mean(plank_pixels)
        threshold = thresh(plank_avg, fixed_thresh)
        mask = img .>= threshold
		masks[:,:,t] = mask
	end
end

function output_images!(stack, masks, overlay, dir, well)
    flat_stack = vec(stack)
    img_min = quantile(flat_stack, 0.0035)
    img_max = quantile(flat_stack, 0.9965)
    adjust_histogram!(stack, LinearStretching(src_minval=img_min, src_maxval=img_max, 
                                              dst_minval=0, dst_maxval=1))
    save("$dir/results_images/$well"*".tif", stack)
    @inbounds for i in CartesianIndices(stack)
        gray_val = RGB{N0f8}(stack[i], stack[i], stack[i])
        overlay[i] = masks[i] ? RGB{N0f8}(1, 0, 0) : gray_val
    end
    save("$dir/results_images/$well"*"mask.tif", overlay)
end

function main()
    dir  = "/mnt/f/Sandhya_Imaging/Time_Lapses/" 
    images = sort([f for f in readdir(dir) if occursin(".ome.tif", f)], lt=natural)
    blockDiameter = 101 
    shift_thresh = 100 
	if isdir("$dir/results_images")
		rm("$dir/results_images"; recursive = true)
	end
	if isdir("$dir/results_data")
		rm("$dir/results_data"; recursive = true)
	end
	mkdir("$dir/results_images")
	mkdir("$dir/results_data")
	output_file = "$dir/results_data/BF_imaging.csv"
	data_matrix = Array{Float64, 2}(undef, ntimepoints, num_wells)
	all_wells = vcat(values(conditions)...)
	conditions = values(conditions)
	@inbounds for (i, wells) in enumerate(conditions)
		@inbounds for j in eachindex(wells)
			images = Array{Gray{N0f16}, 3}(undef, height, width, ntimepoints)
			well = wells[j]
			well_files = sort([f for f in readdir(dir) if occursin(well*"_", f)], 
						 lt=natural)
			read_images!(well, dir, height, width, ntimepoints, images, well_files)
			images = Float64.(images)
			normalized_stack = similar(images)
			registered_stack = similar(images)
			fpMax = maximum(images) 
			fpMin = minimum(images) 
			fpMean = (fpMax - fpMin) / 2.0 + fpMin
			fixed_thresh = fpMean - 0.04
			images, output_stack = stack_preprocess(images, normalized_stack,
													registered_stack, blockDiameter,
													shift_thresh, fpMean, ntimepoints, 
													shift_thresh, sig)
			masks = zeros(Bool, size(images))
			compute_mask!(output_stack, images, masks, sig, 
						  fixed_thresh, ntimepoints)
			output_stack = Gray{N0f8}.(output_stack)
			overlay = zeros(RGB{N0f8}, size(output_stack)...)
			output_images!(output_stack, masks, overlay, dir, well)
			let images = images
				@floop for t in 1:ntimepoints
					@inbounds signal = @views mean((1 .- images[:,:,t]) .* masks[:,:,t])
					@inbounds data_matrix[t, i*j] = signal 
				end
			end
		end # loop over wells for a condition 
	end # loop over conditions 
	df = DataFrame(data_matrix, Symbol.(all_wells))
	df .= ifelse.(isnan.(df), 0, df)
	write(output_file, df)
end
main()
