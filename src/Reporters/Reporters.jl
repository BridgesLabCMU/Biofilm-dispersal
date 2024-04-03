using Images: imfilter, mapwindow, adjust_histogram!, LinearStretching, Gray, N0f16, N0f8,
              Kernel, warp, axes, KernelFactors, RGB 
using ImageMorphology: label_components, component_lengths
using ImageView
using StatsBase: mean, median, quantile
using TiffImages
using FileIO
using NaturalSort: sort, natural
using DataFrames: DataFrame
using CSV: write
using JSON: parsefile
using IntegralArrays: IntegralArray
using CoordinateTransformations: Translation
using IntervalSets: width, leftendpoint, rightendpoint, Interval, ±
using AbstractFFTs
using FFTW
using Compat
using FLoops

thresh(μ, t₀) = (t₀-0.02) + ((t₀+0.06)-(t₀-0.06))*μ

function phase_offset(source::AbstractArray, target::AbstractArray; kwargs...)
    plan = plan_fft(source)
    return phase_offset(plan, plan * source, plan * target; kwargs...)
end

function phase_offset(
    plan,
    source_freq::AbstractMatrix{<:Complex{T}},
    target_freq;
    upsample_factor = 1,
    normalize = false,
) where {T}
    image_product = @. source_freq * conj(target_freq)
    if normalize
        @. image_product /= max(abs(image_product), eps(T))
    end
    if isone(upsample_factor)
        cross_correlation = ifft!(image_product)
    else
        cross_correlation = plan \ image_product
    end
    maxima, maxidx = @compat findmax(abs, cross_correlation)
    shape = size(source_freq)
    midpoints = map(ax -> (first(ax) + last(ax)) / T(2), axes(source_freq))
    idxoffset = map(first, axes(cross_correlation))
    shift = @. T(ifelse(maxidx.I > midpoints, maxidx.I - shape, maxidx.I) - idxoffset)

    isone(upsample_factor) &&
        return (; shift, calculate_stats(maxima, source_freq, target_freq)...)

    shift = @. round(shift * upsample_factor) / T(upsample_factor)
    upsample_region_size = ceil(upsample_factor * T(1.5))
    dftshift = div(upsample_region_size, 2)
    sample_region_offset = @. dftshift - shift * upsample_factor
    cross_correlation = upsampled_dft(
        image_product,
        upsample_region_size,
        upsample_factor,
        sample_region_offset,
    )
    maxima, maxidx = @compat findmax(abs, cross_correlation)
    shift = @. shift + (maxidx.I - dftshift - idxoffset) / T(upsample_factor)

    stats = calculate_stats(maxima, source_freq, target_freq)
    return (; shift, stats...)
end

function upsampled_dft(
    data::AbstractMatrix{T},
    region_size,
    upsample_factor,
    offsets,
) where {T<:Complex}
    shiftrange = 1:region_size
    idxoffset = map(first, axes(data))
    sample_rate = inv(T(upsample_factor))
    freqs = fftfreq(size(data, 2), sample_rate)
    kernel = @. cis(-T(2π) * (shiftrange - offsets[2] - idxoffset[2]) * freqs')

    _data = kernel * data'

    freqs = fftfreq(size(data, 1), sample_rate)
    kernel = @. cis(T(2π) * (shiftrange - offsets[1] - idxoffset[1]) * freqs')
    _data = kernel * _data'
    return _data
end

function calculate_stats(crosscor_maxima, source_freq, target_freq)
    source_amp = mean(abs2, source_freq)
    target_amp = mean(abs2, target_freq)
    error = 1 - abs2(crosscor_maxima) / (source_amp * target_amp)
    phasediff = atan(imag(crosscor_maxima), real(crosscor_maxima))
    return (; error, phasediff)
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

function stack_preprocess(img_stack, normalized_stack, registered_stack,
                        shift_thresh, fpMean, nframes, mxshift, sig, constitutive, blockDiameter)       
    shifts = (0.0, 0.0) 
    @inbounds for t in 1:nframes
        img = img_stack[:,:,constitutive,t]
        img_copy = copy(img)
        img_normalized = normalize_local_contrast(img, img_copy, blockDiameter, fpMean)
        normalized_stack[:,:,t] = imfilter(img_normalized, Kernel.gaussian(sig))
        if t == 1
            registered_stack[:,:,t] = normalized_stack[:,:,t]
        else
            moving = normalized_stack[:,:,t]
            fixed = normalized_stack[:,:,t-1]
            shift, _, _ = phase_offset(fixed, moving, upsample_factor=1)
            if sqrt(shift[1]^2 + shift[2]^2) >= mxshift
                registered_stack[:,:,t] = warp(moving, Translation(shifts[1], shifts[2]), axes(fixed))
                img_stack[:,:,:,t] = warp(img_stack[:,:,:,t], Translation(shifts[1], shifts[2], 0), axes(fixed))
            else
                shift = Tuple([-1*shift[1], -1*shift[2]])
                shift = shift .+ shifts
                shifts = shift
                registered_stack[:,:,t] = warp(moving, Translation(shift[1], shift[2]), axes(fixed))
                img_stack[:,:,:,t] = warp(img_stack[:,:,:,t], Translation(shift[1], shift[2], 0), axes(img_stack[:,:,:,t]))
            end
        end
    end
    processed_stack, crop_indices = crop(registered_stack)
    row_min, row_max, col_min, col_max = crop_indices
    img_stack = img_stack[row_min:row_max, col_min:col_max, :, :]
    return img_stack, processed_stack
end

function compute_mask!(stack, raw_stack, masks, 
                        sig, fixed_thresh, ntimepoints, constitutive)
    @inbounds for t in 1:ntimepoints
        @views img = stack[:,:,t]
        plank_mask = img .< fixed_thresh 
        @views plank = plank_mask .* raw_stack[:,:,constitutive,t] 
        flattened_plank = vec(plank)
        plank_pixels = filter(x -> x != 0, flattened_plank)
        plank_avg = mean(plank_pixels)
        threshold = thresh(plank_avg, fixed_thresh)
        mask = img .>= threshold
		masks[:,:,t] = mask
	end
end

function output_images!(stack, masks, overlay, output_dir)
    imshow(masks)
    flat_stack = vec(stack)
    img_min = quantile(flat_stack, 0.0035)
    img_max = quantile(flat_stack, 0.9965)
    adjust_histogram!(stack, LinearStretching(src_minval=img_min, src_maxval=img_max, 
                                              dst_minval=0, dst_maxval=1))
    TiffImages.save("$output_dir/constitutive.tif", stack)
    @inbounds for i in CartesianIndices(stack)
        gray_val = RGB{N0f8}(stack[i], stack[i], stack[i])
        overlay[i] = masks[i] ? RGB{N0f8}(1, 0, 0) : gray_val
    end
    TiffImages.save("$output_dir/constitutive_mask.tif", overlay)
end

function main()
    dir = "/mnt/f/Sandhya_Imaging/Time_Lapses/Reporters/" 
    files = sort([f for f in readdir(dir, join=true) if occursin(".ome.tif", f) && occursin("check", f)], lt=natural)
    shift_thresh = 100 
    sig = 2 
    constitutive = 2 # index for constitutive channel
    reporter = 1 # index for reporter channel
    blockDiameter = 301
    for file in files
        if isdir("$dir/results_images_"*basename(file)[1:end-8])
            rm("$dir/results_images_"*basename(file)[1:end-8]; recursive = true)
        end
        if isdir("$dir/results_data_"*basename(file)[1:end-8])
            rm("$dir/results_data_"*basename(file)[1:end-8]; recursive = true)
        end
        output_img_dir = "$dir/results_images_"*basename(file)[1:end-8]
        output_data_dir = "$dir/results_data_"*basename(file)[1:end-8]
        mkdir(output_img_dir)
        mkdir(output_data_dir)
        output_file = "$output_data_dir/data.jld2"
        height, width, channels, nslices, ntimepoints = size(FileIO.load(file))
        @views images = Float64.(sum(reshape(TiffImages.load(file), 
                                             (height, width, nslices, channels, ntimepoints)), dims=3)[:,:,1,:,:])
        data_matrix = Array{Float64, 1}(undef, ntimepoints)
        normalized_stack = images[:,:,constitutive,:]
        registered_stack = copy(normalized_stack)
        fpMax = maximum(normalized_stack) 
        fpMin = minimum(normalized_stack) 
        fpMean = (fpMax - fpMin) / 2.0 + fpMin
        fixed_thresh = fpMean + 0.02
        images, output_stack = stack_preprocess(images, normalized_stack,
                                                registered_stack, 
                                                shift_thresh, fpMean, ntimepoints, 
                                                shift_thresh, sig, constitutive, blockDiameter)
        @views masks = zeros(Bool, size(images[:,:,constitutive,:]))
        compute_mask!(output_stack, images, masks, sig, 
                      fixed_thresh, ntimepoints, constitutive)
        output_stack = Gray{N0f8}.(output_stack)
        overlay = zeros(RGB{N0f8}, size(output_stack)...)
        output_images!(output_stack, masks, overlay, output_img_dir)
        let images = images
            @floop for t in 1:ntimepoints
                @inbounds signal = @views mean((images[:,:,reporter,t] ./ images[:,:,constitutive,t]) .* masks[:,:,t])
                @inbounds data_matrix[t] = signal 
            end
        end
        data_matrix .= ifelse.(isnan.(data_matrix), 0, data_matrix)
        FileIO.save(output_file, Dict("data" => data_matrix))
    end
end
main()
