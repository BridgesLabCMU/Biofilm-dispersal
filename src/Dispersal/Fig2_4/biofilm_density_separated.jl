using Makie
using GLMakie
using CairoMakie
CairoMakie.activate!(type = "pdf")
set_theme!(fonts=(; bold="TeX Gyre Heros Makie"))
using LaTeXStrings
using TiffImages
using NaturalSort: sort, natural
using StatsBase
using ColorTypes: Gray, N0f16
using ImageMorphology: label_components, component_lengths, component_centroids
using Images: imresize, distance_transform, feature_transform
using HistogramThresholding: find_threshold, Otsu
using FLoops
using DelimitedFiles

round_up(x, multiple) = ceil(x / multiple) * multiple

function radial_averaging(data_files, bin_interval)
    first_image = TiffImages.load(data_files[1]) .> 0
    labels = label_components(first_image)
    volumes = component_lengths(labels)
    centers = component_centroids(labels)
    center = centers[argmax(volumes[1:end])]
    center = (round(Int, center[1]), round(Int, center[2]), 1)
    center_distance = zeros(Bool, size(labels))
    center_distance[center[1], center[2], center[3]] = true
    center_distance = distance_transform(feature_transform(center_distance))
    max_distance = round_up(maximum(center_distance), bin_interval)
    bins = 0:bin_interval:max_distance
    nbins = length(bins)
    data_distribution = zeros(nbins-1)
    random_distribution = zeros(nbins-1)
    data_configuration = TiffImages.load(data_files[end]) .> 0
    N_voxels = sample(findall(!iszero, first_image), sum(first_image) - sum(data_configuration), replace=false)
    random_configuration = copy(first_image)
    random_configuration[N_voxels] .= 0
    for i in 1:length(data_distribution) 
        data_distribution[i] = mean(data_configuration[findall(x -> bins[i] <= x <= bins[i+1], center_distance)])
        random_distribution[i] = mean(random_configuration[findall(x -> bins[i] <= x <= bins[i+1], center_distance)])
    end
    return data_distribution, random_distribution 
end

function main()
    plot_ylabel = "Density (a.u.)"
    plot_xlabel = "Distance from center (µm)"
    master_directory = "/mnt/h/Dispersal"
    image_folders = filter(isdir, readdir(master_directory, join=true))
    image_folders = [f for f in image_folders if occursin("WT", f)]
    filter!(folder->folder≠master_directory*"/Plots", image_folders)
    plots_folder = "/mnt/h/Dispersal/Plots"
    n = 5 
    xtick_interval = n/(0.065*4*8)

    for (i,images_folder) in enumerate(image_folders)
        plot_filename = "final_timepoint_density_"*string(i) 
        files = sort([f for f in readdir(images_folder, join=true) if occursin("downsampled_mask", f)], 
                                 lt=natural)
        data_distribution, random_distribution = radial_averaging(files, 8)
        xs = 0:xtick_interval:length(data_distribution)-1
        fig = Figure(size=(7*72, 3*72))
        ax = Axis(fig[1, 1])
        lines!(ax, 0:length(data_distribution)-1, data_distribution, label="Data", color=:black, linewidth=2)
        lines!(ax, 0:length(data_distribution)-1, random_distribution, label="Random model",color=:black, linestyle=:dash,linewidth=2)
        writedlm("$(plots_folder)/WT_density_$(string(i)).csv", data_distribution, ",")
        writedlm("$(plots_folder)/WT_random_density_$(string(i)).csv", random_distribution, ",")
        ax.xticks = xs
        ax.xtickformat=values->string.([round(Int, v*n/xtick_interval) for v in values])
        ax.xlabel = plot_xlabel
        ax.ylabel = plot_ylabel
        ax.rightspinevisible = false
        ax.topspinevisible = false
        ax.xgridvisible = false
        ax.ygridvisible = false
        fig[1,2] = Legend(fig, ax, framevisible=false, labelsize=12, rowgap=0)
        save("$plots_folder/$plot_filename"*".pdf", fig)
    end
end

main()
