function strel_circle(radius)
    # Circular structuring element
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

function classify_pixels(curr_mask, prev_mask; thresh=0.016)
    # Identify pixels that are likely to belong to the biofilm based on distance from the previous mask
    distance = distance_transform(feature_transform(prev_mask)) .* curr_mask 
    normalized_distance = distance ./ maximum(vec(distance))
    refined_pixels = (normalized_distance .< thresh) .* curr_mask 
    return refined_pixels
end

function mask_prep!(prev_mask, curr_mask, prev_img, curr_img, t, intensity_thresholds, height, width, slices)
    # Generate mask using basic intensity threshold
    @floop for i in 1:slices
        @views curr_mask[:,:,i] = curr_img[:,:,i] .> intensity_thresholds[i]
    end
    return nothing 
end

function isolate_biofilm!(curr_mask, prev_mask, closed_prev, slices, t, timepoint_thresh)
    # Remove planktonic cells at the current timepoint using morphological reconstruction on masks of current
    # and previous cells
	if t < timepoint_thresh || timepoint_thresh < 2
		union = prev_mask .| curr_mask 
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
	else
        # Once past a specified timepoint, just refine the current mask using the previous mask and
        # a basic classification of pixels likely to belong to the biofilm (i.e., those that are close
        # to the previous mask)
        curr_mask .*= (closed_prev .* classify_pixels(curr_mask, prev_mask))
    end
end
