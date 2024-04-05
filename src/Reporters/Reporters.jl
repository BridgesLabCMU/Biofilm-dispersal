using Images: imfilter, mapwindow, adjust_histogram!, LinearStretching, Gray, N0f16, N0f8,
              Kernel, warp, axes, KernelFactors, RGB 
using ImageView
using StatsBase: mean, median, quantile
using TiffImages
using FileIO
using NaturalSort: sort, natural
using IntegralArrays: IntegralArray
using CoordinateTransformations: Translation
using IntervalSets: width, leftendpoint, rightendpoint, Interval, ±
using AbstractFFTs
using FFTW
using Compat
using FLoops
using PythonCall

morph = pyimport("skimage.morphology")

#thresh(μ, t₀) = ((t₀+0.06)-(t₀-0.06))*μ

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

function stack_preprocess(img_stack, normalized_stack, registered_stack,
                        shift_thresh, fpMean, nframes, mxshift, constitutive, blockDiameter)       
    shifts = (0.0, 0.0) 
    @inbounds for t in 1:nframes
        img = img_stack[:,:,constitutive,t]
        img_normalized = img .- mean(img)
        img_normalized .+= fpMean
        normalized_stack[:,:,t] = img_normalized
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

function compute_mask!(stack, masks, fixed_thresh, ntimepoints, constitutive)
    @inbounds for t in 1:ntimepoints
        @views img = mapwindow(median, stack[:,:,t], (5,5))
        flat_img = vec(img)
        img_min = quantile(flat_img, 0.0035)
        img_max = quantile(flat_img, 0.9965)
        adjust_histogram!(img, LinearStretching(src_minval=img_min, src_maxval=img_max, 
                                                  dst_minval=0, dst_maxval=1))
        mask = img .> fixed_thresh 
        mask = pyconvert(Array, morph.isotropic_opening(mask, radius=5))
		masks[:,:,t] = mask
	end
end

function output_images!(stack, masks, overlay, output_dir, replicate)
    imshow(masks)
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
    dir = "/mnt/f/Sandhya_Imaging/Time_Lapses/Reporters/Ptac/" # Folder containing the images for strain X 
    files = sort([f for f in readdir(dir, join=true) if occursin(".ome.tif", f)], lt=natural)
    shift_thresh = 100 
    sig = 2 
    constitutive = 2 # index for constitutive channel
    reporter = 1 # index for reporter channel
    blockDiameter = 301
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
        normalized_stack = images[:,:,constitutive,:]
        registered_stack = copy(normalized_stack)
        fpMax = maximum(normalized_stack) 
        fpMin = minimum(normalized_stack) 
        fpMean = (fpMax - fpMin) / 2.0 + fpMin
        images, output_stack = stack_preprocess(images, normalized_stack,
                                                registered_stack, 
                                                shift_thresh, fpMean, nframes, 
                                                shift_thresh, constitutive, blockDiameter)
        @views first_slice = output_stack[:,:,1]
        output_stack = output_stack .- first_slice[:,:,:]
        output_stack[output_stack .< 0] .= 0
        @views masks = zeros(Bool, size(output_stack))
        fpMax = maximum(output_stack) 
        fpMin = minimum(output_stack) 
        fpMean = (fpMax - fpMin) / 2.0 + fpMin
        fixed_thresh = fpMean + 0.01
        @show fixed_thresh
        compute_mask!(output_stack, masks, 
                      fixed_thresh, nframes, constitutive)
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
