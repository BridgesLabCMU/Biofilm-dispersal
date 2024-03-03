module RemovePlanktonicCells

using Revise
using ImageMorphology: mreconstruct, dilate, closing, opening, 
                            strel, label_components, component_lengths
using Images: distance_transform, feature_transform
using ColorTypes: Gray, N0f16
using TiffImages: save, load
using StatsBase: mean, std
using FLoops

export isolate_biofilm!, strel_circle, mask_prep_prethresh!, mask_prep_postthresh 

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

function classify_pixels(curr_mask, curr_cv, prev_mask; w1=0.9, w2=0.1,
                            thresh=0.98)
    distance = distance_transform(feature_transform(.!prev_mask)) .* curr_mask 
    curr_cv .*= curr_mask 
    normalized_distance = distance ./ maximum(vec(distance))
    normalized_pixel_values = curr_cv ./ maximum(vec(curr_cv))
    probability = w1 .* (1 .- normalized_distance) + w2 .* (1 .- normalized_pixel_values)
    refined_pixels = (probability .> thresh) .* curr_mask 
    return refined_pixels
end

function mask_prep_prethresh!(prev_mask, curr_mask, prev_img, curr_img, t, intensity_thresholds, height, width, slices)
    if t == 2
        @floop for i in 1:slices
            prev_mask[:,:,i] = closing(prev_img[:,:,i] .> intensity_thresholds[i], 
                                              strel_circle(6))
            @views curr_mask[:,:,i] = curr_img[:,:,i] .> intensity_thresholds[i]
        end
    else
        @floop for i in 1:slices
            @views curr_mask[:,:,i] = curr_img[:,:,i] .> intensity_thresholds[i]
        end
    end
    return nothing 
end

function mask_prep_postthresh(curr_mask, curr_img, t, files, dir, intensity_thresholds, height, width, slices, frames, half_window)
    first_index = max(1, t-half_window)
    last_index = min(frames, t+half_window)
    window = Array{N0f16, 4}(undef, height, width, slices, last_index-first_index+1) 
    window_mask = zeros(Bool, (height, width, slices, last_index-first_index+1))
   
    @floop for i in 1:slices
        @views curr_mask[:,:,i] = curr_img[:,:,i] .> intensity_thresholds[i]
    end
    
    @inbounds for τ in first_index:last_index
        @views file = files[τ]
        window[:,:,:,τ-first_index+1] = load("$dir/$file")
    end
    
    @inbounds for τ in size(window, 4)
        @inbounds for i in 1:slices
            window_mask[:,:,i,τ] = @views window[:,:,i,τ] .> intensity_thresholds[i]
        end
    end

    window = nothing
    mask_float = Float32.(window_mask)
    window_mask = nothing
    mask_mean = mean(mask_float, dims=4) 
    mask_std = std(mask_float, dims=4) 
    @views curr_cv = @. ifelse(iszero(mask_mean), 0.0, mask_std/mask_mean)[:,:,:,1]
    mask_float = nothing
    mask_mean = nothing
    mask_std = nothing
    return curr_img, curr_mask, curr_cv
end

function isolate_biofilm!(curr_mask, prev_mask, curr_cv, slices)
	if curr_cv == nothing 
		union = prev_mask .| curr_mask 
	else
		union = prev_mask .| classify_pixels(curr_mask, curr_cv, prev_mask) 
	end
	curr_recon = mreconstruct(dilate, prev_mask, union)
	curr_recon .*= curr_mask 
	@inbounds for i in 1:slices
		closed = closing(curr_recon[:,:,i], strel_circle(6))
		opened = opening(closed, strel_circle(5))
		@views absent = label_components(opened .== false .& curr_mask[:,:,i])
		areas = component_lengths(absent)
        @views absent[absent .∈ (findall(areas[1:end] .> 300),)] .= 0
        @views absent[absent .∈ (findall(areas[1:end] .< 300),)] .= 1
        union = absent .| opened 
		curr_recon[:,:,i] = mreconstruct(dilate, opened, union)
	end
	curr_mask .*= curr_recon
end

end # module
