module RemovePlanktonicCells

using Revise
using ImageMorphology: mreconstruct, dilate, closing, opening, 
                            strel, label_components, component_lengths
using Images: distance_transform, feature_transform
using ColorTypes: Gray, N0f16

export isolate_biofilm!, strel_circle 

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

function classify_pixels(curr_mask, curr_cv, prev_mask; w1=0.8, w2=0.2,
                            thresh=0.98)
    distance = distance_transform(feature_transform(.!prev_mask)) .* curr_mask 
    curr_cv .*= curr_mask 
    normalized_distance = Float32.(distance) ./ maximum(vec(distance))
    normalized_pixel_values = curr_cv ./ maximum(vec(curr_cv))
    probability = @. w1 * (1 - normalized_distance) + w2 * (1 - normalized_pixel_values)
    refined_pixels = @. (probability > thresh) * curr_mask 
    return refined_pixels
end

function isolate_biofilm!(curr_img, curr_mask, prev_mask, curr_cv, slices)
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
	curr_img .*= curr_mask 
end

end # module
