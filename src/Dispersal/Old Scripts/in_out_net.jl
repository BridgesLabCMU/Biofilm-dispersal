using Plots
using LaTeXStrings
using TiffImages: load, save
using NaturalSort: sort, natural
using StatsBase
using ColorTypes: Gray, N0f16
using ImageMorphology: label_components, component_lengths, component_centroids
using Images: imresize, distance_transform, feature_transform
using HistogramThresholding: find_threshold, Otsu
using FLoops
using DelimitedFiles

pgfplotsx()

round_up(x, multiple) = ceil(x / multiple) * multiple

function N_smallest(M, center_distance, biofilm_configuration, n)
    M[center_distance .== 0] .= biofilm_configuration[center_distance .== 0]  
    v = vec(M)
    nonzero_indices = findall(x -> x != 0, v)
    n = min(n, length(nonzero_indices))
    if n == 0
        return CartesianIndex[]
    end
    perm = partialsortperm(v[nonzero_indices], 1:n)
    indices = CartesianIndices(M)[nonzero_indices[perm]]
    return indices
end

function radial_averaging(files, images_folder, first_index, end_index, bin_interval, budget)
    first_image = load("$images_folder/$(files[first_index])")
    labels = label_components(first_image)
    volumes = component_lengths(labels)
    centers = component_centroids(labels)
    center = centers[argmax(volumes[1:end])]
    center = (round(Int, center[1]), round(Int, center[2]), 10)
    center_distance = zeros(Bool, size(labels))
    center_distance[center[1], center[2], center[3]] = true
    center_distance = distance_transform(feature_transform(center_distance))
    max_distance = round_up(maximum(center_distance), bin_interval)
    bins = 0:bin_interval:max_distance
    nbins = length(bins)
    ntimepoints = end_index - first_index
    data_matrix = zeros(nbins-1, ntimepoints)
    biofilm_configuration = load("$images_folder/$(files[first_index])")
    dispersal_config = zeros(Bool, size(biofilm_configuration))
    for t in 1:ntimepoints
        dispersal_config .= 0 
        N_voxels = N_smallest(biofilm_configuration .* center_distance, center_distance, 
                              biofilm_configuration, budget[t])
        if N_voxels == CartesianIndex[]
            for i in 1:size(data_matrix, 1) 
                data_matrix[i, t] = 0 
            end
        else
            biofilm_configuration[N_voxels] .= 0 
            dispersal_config[N_voxels] .= 1 
            for i in 1:size(data_matrix, 1) 
                data_matrix[i, t] = -1*mean(dispersal_config[findall(x -> bins[i] <= x <= bins[i+1], center_distance)])
            end
        end
    end
    return data_matrix
end

function main()
    push!(PGFPlotsX.CUSTOM_PREAMBLE, "\\usepackage{sfmath}\n\\renewcommand{\\familydefault}{\\sfdefault}")
    default(titlefont = 20, legendfontsize = 15, 
            guidefont = (20, :black), colorbar_tickfontsize=15, colorbar_titlefontsize=20, tickfont = (15, :black), 
            guide = L"x", linewidth=2, grid=false, formatter=:plain)
    plot_size = (350,300)
    plot_xlabel = "Time (h)"
    plot_ylabel = L"Distance from center ($\mu$m)"
    plot_title = "Inside-out model"

    master_directory = "/mnt/h/Dispersal"
    image_folders = filter(isdir, readdir(master_directory, join=true))
    filter!(folder->folder≠master_directory*"/Plots", image_folders)
    plots_folder = "/mnt/h/Dispersal/Plots"

    for images_folder in image_folders
        plot_filename = basename(images_folder)*"_in_out_net" 
        if isfile("$plots_folder/$plot_filename"*".svg")
            continue
        end
        files = sort([f for f in readdir(images_folder) if occursin("mask_isotropic", f)], 
                                 lt=natural)
        ntimepoints = length(files)
        first_image = load("$images_folder/$(files[1])"; lazyio=true)
        height, width, slices = size(first_image)
        first_image = nothing
        net = readdlm(plots_folder*"/"*basename(images_folder)*".csv", ',', Int)[1:end,1]
        first_index = argmax(net)
        end_index = min(first_index + 45, ntimepoints)
        budget = -1 .* diff(net[first_index:end_index])
        budget[budget .< 0] .= 0
        data_matrix = radial_averaging(files, images_folder, first_index, end_index, 30, budget)
        replace!(data_matrix, Inf=>NaN)
        replace!(data_matrix, NaN=>0.0)
        writedlm("$(plots_folder)/$(plot_filename).csv", data_matrix, ",")
        n = 10
        ytick_interval = n/0.065/30
        xs = 0:6:size(data_matrix, 2)-1
        ys = 0:ytick_interval:size(data_matrix, 1)-1
        c = cgrad(:Purples, rev=true)
        plt = heatmap(data_matrix, xticks=xs, yticks=ys, color=c, clim=(-0.1, 0),
                      colorbar_title="Density change (a.u.)", colorbar_ticks=(-0.1:0.02:0), xformatter=xi -> xi*1/6, 
                      yformatter=yi -> yi/ytick_interval*n, size=plot_size)
        xlabel!(plt, plot_xlabel)
        ylabel!(plt, plot_ylabel)
        title!(plt, plot_title)
        savefig("$plots_folder/$plot_filename"*".svg")
    end
end

main()
