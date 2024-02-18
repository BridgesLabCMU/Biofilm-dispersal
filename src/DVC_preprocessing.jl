using Images: warp, channelview, Gray, N0f16, mapwindow, distance_transform,
              feature_transform
using ImageTransformations: imresize
using Interpolations: Linear
using StatsBase: mean, std
using TiffImages: load, save, IMAGEDESCRIPTION, ifds
using ImageMorphology: tophat!, mreconstruct, dilate, closing, opening, 
                            strel, label_components, component_lengths
using CoordinateTransformations: Translation
using HistogramThresholding: find_threshold, Otsu
using FLoops
using Revise
using AbstractFFTs
using Compat
using FFTW
using Statistics
using TensorOperations
using ProgressMeter

function rolling_mean(X::AbstractArray{T, 4}, window_size::Int) where {T<:Number}
    if iseven(window_size)
        error("Window size must be odd")
    end
	half_window = div(window_size, 2) 
	ntimepoints = size(X, 4)
	Y = similar(X)
	@inbounds for n in 1:ntimepoints
		lo = max(1, n-half_window)
		hi = min(ntimepoints, n+half_window)
		Y[:,:,:,n] = mean(X[:,:,:,lo:hi], dims=4)
	end
	return Y
end

function rolling_std(X::AbstractArray{T, 4}, window_size::Int) where {T<:Number}
    if iseven(window_size)
        error("Window size must be odd")
    end
	half_window = div(window_size, 2) 
	ntimepoints = size(X, 4)
	Y = similar(X)
	@inbounds for n in 1:ntimepoints
		lo = max(1, n-half_window)
		hi = min(ntimepoints, n+half_window)
		Y[:,:,:,n] = std(X[:,:,:,lo:hi], dims=4)
	end
	return Y
end

function strel_circle(radius)
    diameter = 2 * radius + 1
    mask = zeros(Bool, (diameter, diameter))
    @inbounds for i in 1:diameter
        @inbounds for j in 1:diameter
            if (i-radius-1)^2 + (j-radius-1)^2 <= radius^2
                mask[i,j] = true
            end
        end
    end
    return strel(Bool, mask)
end


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
    cropped_stack = img_stack[i1:i2, j1:j2, z1:z2, :]
    return cropped_stack
end

function register!(img_stack, registered_stack, nframes, center)       
    reference = img_stack[:,:,:,center]
    @floop for t in 1:nframes
        moving = img_stack[:,:,:,t]
        shift, _, _ = phase_offset(reference, moving, upsample_factor=10)
        shift = Translation(-1*shift[1], -1*shift[2], -1*shift[3])
        registered_stack[:,:,:,t] = warp(img_stack[:,:,:,t], 
                                         shift, axes(registered_stack[:,:,:,t]), 1)
    end
    return nothing
end

function timepoint_mask!(timeseries_mask, first_timepoints, nframes)
    @inbounds for t in 1:nframes 
        @views mask = timeseries_mask[:,:,:,t] .& (first_timepoints .== -1)
        first_timepoints[mask] .= t
    end
    return nothing
end

function fill_mask!(mask, timeseries, slices, cell_threshold)
    @floop for i in 1:slices 
        @views slice_ = timeseries[:, :, i, :]
        otsu_thresh = find_threshold(slice_, Otsu())
        if otsu_thresh/3 < cell_threshold 
            mask[:, :, i, :] = slice_ .> cell_threshold 
        else
            mask[:, :, i, :] = slice_ .> otsu_thresh/3
        end
    end
    return nothing 
end

function classify_pixels(curr_mask, std_img, prev_mask; w1=0.8, w2=0.2,
                            thresh=0.98)
    distance = distance_transform(feature_transform(~prev_mask .* curr_mask)) 
    std_image .*= curr_mask 
    normalized_distance = reinterpret(Float64, distance) ./ max.(distance)
    normalized_pixel_values = std_image ./ max.(std_image)
    probability = @. w1 * (1 - normalized_distance) + w2 * (1 - normalized_pixel_values)
    refined_pixels = @. (probability > thresh) * curr_mask 
    return refined_pixels
end

function isolate_biofilm!(timeseries, height, width, slices, frames, cell_threshold, dir)
    mask = zeros(Bool, (height, width, slices, frames))
    fill_mask!(mask, timeseries, slices, cell_threshold)
    mask_float = Float64.(mask)
    mask_mean = rolling_mean(mask_float, 5) 
    mask_std = rolling_std(mask_float, 5) 
    std_img = @. ifelse(iszero(mask_mean), mask_std/mask_mean, 0.0)
    first_timepoints = fill(-1, (height, width, slices))
    timepoint_mask!(mask, first_timepoints, frames)
    timepoint_thresh = find_threshold(first_timepoints, Otsu())
    prev_mask = zeros(Bool, (height, width, slices)) 
    @floop for i in 1:slices
        @views prev_mask[:,:,i] = closing(mask[:,:,i,1], strel_circle(6))
    end
    timeseries[:,:,:,1] .*= prev_mask
    @showprogress dt=1 desc="Removing planktonic cells..." for t in 2:frames
        if t < timepoint_thresh
            @views union = prev_mask .| mask[:,:,:,t]
        else
            @views union = prev_mask .| classify_pixels(mask[:,:,:,t], std_img[:,:,:,t], prev_mask) 
        end
        curr_recon = mreconstruct(dilate, prev_mask, union)
        @views curr_recon .*= mask[:,:,:,t]
        @inbounds for i in 1:slices
            closed = closing(curr_recon[:,:,i,t], strel_circle(6))
            opened = opening(closed, strel_circle(5))
            @views absent = label_components(opened .== false .& mask[:,:,i,t])
            areas = component_lengths(absent)
            @views union = absent[absent .∈ (findall(areas[1:end] .> 300),)] .| opened 
            curr_recon[:,:,i,t] = mreconstruct(dilate, opened, union)
        end
        curr_mask .*= curr_recon
        timeseries[:,:,:,t] .*= curr_mask 
        prev_mask = curr_mask 
    end
end

function isotropic_voxel(zstack, aspect_ratio)
    return imresize(zstack, ratio=(1,1,aspect_ratio), method=Linear())
end

function write_zstacks!(timeseries, frames, dir, aspect_ratio)
    @floop for t in 1:frames
        #isotropic = isotropic_voxel(timeseries[:,:,:,t], aspect_ratio)
        isotropic = timeseries[:,:,:,t]
        save("$dir/stack_$(t).tif", isotropic)
    end
end

function main()
    cell_threshold = 3.0e-5 # Depends on imaging setup!!
    file = "cheY_replicate1_noback.tif" 
    if isdir(file[1:end-10]*"processed")
        rm(file[1:end-10]*"processed"; recursive = true)
	end
	mkdir(file[1:end-10]*"processed") 
    dir = file[1:end-10]*"processed"
    timeseries = load(file)
    height, width, zt = size(timeseries)
    ImageJ_metadata = first(ifds(timeseries))[IMAGEDESCRIPTION].data
    xy_res = 0.065
    z_spacing = findfirst("spacing=", ImageJ_metadata)
    nloop = findfirst("loop=", ImageJ_metadata)
    z_res = tryparse(Float64, ImageJ_metadata[(z_spacing[8]+1):(nloop[1]-1)])
    aspect_ratio = z_res / xy_res
    slice_loc = findfirst("slices=", ImageJ_metadata)
    frames_loc = findfirst("frames=", ImageJ_metadata)
    hyperstack_loc = findfirst("hyperstack=", ImageJ_metadata)
    slices = tryparse(Int, ImageJ_metadata[(slice_loc[7]+1):(frames_loc[1]-1)])
    frames = tryparse(Int, ImageJ_metadata[(frames_loc[7]+1):(hyperstack_loc[1]-1)])
    timeseries = channelview(reshape(timeseries, (height, width, slices, frames)))
    #noback = tophat!(timeseries, noback; dims=(1,2), r=15)
    registered = similar(timeseries)
    center = frames ÷ 2
    register!(timeseries, registered, frames, center)
    registered = reinterpret(Gray{N0f16}, crop(registered))
    height, width, slices, frames = size(registered)
    isolate_biofilm!(registered, height, width, slices, frames, cell_threshold, dir)
    write_zstacks!(registered, frames, dir, aspect_ratio)
end

main()
