module Preprocessing

using Images: warp 
using TiffImages: save
using ImageMorphology: tophat
using CoordinateTransformations: Translation
using HistogramThresholding: find_threshold, Otsu
using FLoops
using Revise
using AbstractFFTs
using Compat
using FFTW
using Statistics
using TensorOperations

export write_images!, register!, crop, mask_thresholds!, timepoint_threshold

function phase_offset(source, target; kwargs...)
    plan = plan_fft(source)
    return phase_offset(plan, plan * source, plan * target; kwargs...)
end

function phase_offset(plan, source_freq::AbstractArray{<:Complex{T}},
                        target_freq; upsample_factor = 1, 
                        normalize = false) where {T}
    
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

function contract_tensors(kernel, data)
    data = conj.(data)
    if ndims(data) == 2
        @tensor begin
            _data[i,j] := kernel[i,:] * data[j,:]
        end
    else
        @tensor begin
            _data[i,j,k] := kernel[i,:] * data[j,k,:]
        end
    end
    return _data 
end

function upsampled_dft(data::AbstractArray{T}, region_size, upsample_factor,
                            offsets) where {T<:Complex}
    
    shiftrange = 1:region_size
    sample_rate = inv(T(upsample_factor))
    idxoffsets = map(first, axes(data))
    shape = size(data)
    offsets = offsets
    _data = copy(data)

    for (k,(dim, offset, idxoffset)) in enumerate(zip(reverse(shape), 
                                                      reverse(offsets), 
                                                      reverse(idxoffsets)))
        if iseven(k)
            sign = 1
        else
            sign = -1
        end
        freqs = fftfreq(dim, sample_rate)
        kernel = @. cis(sign*T(2π) * (shiftrange - offset - idxoffset) * freqs')
        _data = contract_tensors(kernel, _data) 
    end
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
    @views mask = .!any(isone.(img_stack), dims=4)[:,:,:,1]
    @views mask_i = any(mask, dims=(2,3))[:,1,1]
    @views mask_j = any(mask, dims=(1,3))[1,:,1]
    @views mask_z = any(mask, dims=(1,2))[1,1,:]
    i1 = findfirst(mask_i)
    i2 = findlast(mask_i)
    j1 = findfirst(mask_j)
    j2 = findlast(mask_j)
    z1 = findfirst(mask_z)
    z2 = findlast(mask_z)
    @views cropped_stack = img_stack[i1:i2, j1:j2, z1:z2, :]
    return cropped_stack
end

function register!(img_stack, registered_stack, nframes, center)       
    @views reference = img_stack[:,:,:,center]
    @floop for t in 1:nframes
        @views moving = img_stack[:,:,:,t]
        shift, _, _ = phase_offset(reference, moving, upsample_factor=10)
        shift = Translation(-1*shift[1], -1*shift[2], -1*shift[3])
        registered_stack[:,:,:,t] = warp(moving, shift, axes(moving), 1)
    end
    return nothing
end

function timepoint_mask!(timeseries, first_timepoints, nframes, cell_threshold)
    @inbounds for t in 1:nframes 
        @views mask = (timeseries[:,:,:,t] .> cell_threshold) .& (first_timepoints .== -1)
        first_timepoints[mask] .= t
    end
    return nothing
end

function timepoint_threshold(timeseries, height, width, slices, frames, cell_threshold)
    first_timepoints = fill(-1, (height, width, slices))
    timepoint_mask!(timeseries, first_timepoints, frames, cell_threshold)
    timepoint_thresh = find_threshold(first_timepoints, Otsu())
    return timepoint_thresh
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


function write_images!(timeseries, frames, dir, aspect_ratio)
    @floop for t in 1:frames
        @views isotropic = timeseries[:,:,:,t]
        save("$dir/stack_$(t).tif", isotropic)
    end
end

end # module
