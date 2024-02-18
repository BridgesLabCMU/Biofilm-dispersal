using Plots
using TiffImages: load, save
using ImageSegmentation: feature_transform, distance_transform
using NaturalSort: sort, natural
using ColorTypes: Gray, N0f16
using ImageMorphology: label_components, component_lengths, component_centroids

pgfplotsx()

round_up(x, multiple) = ceil(x / multiple) * multiple

function radial_averaging(binary_timeseries, mask_diff, bin_interval, ntimepoints)
    labels = label_components(binary_timeseries[:,:,:,1])
    volumes = component_lengths(labels)
    centers = component_centroids(labels)
    center = centers[argmax(volumes[1:end])]
    center_distance = distance_transform(feature_transform(center))
    max_distance = round_up(maximum(center_distance), bin_interval)
    bins = 0:bin_interval:max_distance
    nbins = length(bins)
    data_matrix = zeros(nbins, ntimepoints)
    for t in 1:ntimepoints
        for i in eachindex(bins)
            data_matrix[i, t] = mean(mask_diff[center_distance .∈ (bins[i],), t])
        end
    end
    return data_matrix
end

function read_images!(directory, ntimepoints, arr, files)
    @inbounds for t in 1:ntimepoints
        @views file = files[t]
        arr[:, :, :, t] = load("$directory/$file")
    end
    return nothing 
end

function main()
    push!(PGFPlotsX.CUSTOM_PREAMBLE,
    """
    \\usepackage[scaled]{helvet}
    \\renewcommand\\familydefault{\\sfdefault} 
    \\usepackage[T1]{fontenc}
    \\usepackage{helvet, sansmath}
    \\sansmath
    """)
    default(titlefont = 20, legendfontsize = 15, 
            guidefont = (20, :black), colorbar_tickfontsize=15, colorbar_titlefontsize=15, tickfont = (15, :black), 
            guide = L"x", linewidth=2, grid=false, formatter=:plain)

    images_folder = "/Users/jojo/Downloads/FAP_Julia/cheY_replicate1_processed"
    plots_folder = "/Users/jojo/Downloads/FAP_Julia/Plots"
    plot_filename = "cheY_replicate1"
    plot_xlabel = "Time (h)"
    plot_ylabel = L"Distance from center ($\mu$m)"
    files = sort([f for f in readdir(images_folder) if occursin("mask", f)], 
                             lt=natural)
    ntimepoints = length(files)
    dummy_image = load("$images_folder/$(files[1])"; lazyio=true)
    height, width, slices = size(dummy_image)
    timeseries = zeros(Bool, height, width, slices, ntimepoints)
    read_images!(images_folder, ntimepoints, timeseries, files)
    
    mask_diff = diff(reinterpret(Int, timeseries), dims=4) # ntimepoints-1
    net_change = sum(mask_diff, dims=(1:3)) # ntimepoints-1
    dispersal_index = findfirst(x -> x < 0, net_change)
    mask_diff[mask_diff .> 0] .= 0 
    mask_diff = mask_diff[:,:,:,dispersal_index:end]

    data_matrix = radial_averaging(timeseries[:,:,:,dispersal_index:end-1], 
                                   mask_diff, 30, ntimepoints)
	plt = heatmap(plot_xticks, reverse(plot_yticks), data_matrix, 
                  yflip=true, color=:coolwarm,
                  colorbar_title=clab_title, size=plot_size)
    xlabel!(plt, plot_xlabel)
    ylabel!(plt, plot_ylabel)
    title!(plt, plot_title)
    savefig("$plots_folder/$plot_filename"*".svg")
end

main()
