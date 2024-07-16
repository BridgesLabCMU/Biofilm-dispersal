using Makie
using GLMakie
using CairoMakie
CairoMakie.activate!(type = "pdf")
set_theme!(fonts=(; bold="TeX Gyre Heros Makie"))
using SwarmMakie
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
using Trapz

function replace_nan(arr)
    result = copy(arr)
    n = length(arr)
    left = 1
    while left <= n && isnan(arr[left])
        left += 1
    end
    right = n
    while right >= 1 && isnan(arr[right])
        right -= 1
    end
    for i in 1:n
        if isnan(arr[i])
            if i < left
                result[i] = arr[left]
            elseif i > right
                result[i] = arr[right]
            else
                left_diff = i - left
                right_diff = right - i
                if left_diff <= right_diff
                    result[i] = arr[left]
                else
                    result[i] = arr[right]
                end
            end
        end
    end
    return result
end

function centroid(data)
    total_sum = sum(data)
    positions = 1:length(data)
    weighted_sum = sum(data[i] * positions[i] for i in eachindex(data))
    return weighted_sum / total_sum
end

function main()
    plots_folder = "/mnt/h/Dispersal/Plots"
    plot_filenames = ["cheY_centroids", "lapG_centroids", "rbmB_centroids"]
    for plot_filename in plot_filenames
        data_all = []
        in_out_all = []
        out_in_all = []
        random_all = []
        for i in 1:5
            data = readdlm(plots_folder*"/"*plot_filename[1:4]*"_replicate"*string(i)*"_processed_data_intensity.csv", ',')[1:end,1:end,1]
            in_out = readdlm(plots_folder*"/"*plot_filename[1:4]*"_replicate"*string(i)*"_processed_in_out_net_downsampled.csv", ',')[1:end,1:end,1]
            out_in = readdlm(plots_folder*"/"*plot_filename[1:4]*"_replicate"*string(i)*"_processed_out_in_net_downsampled.csv", ',')[1:end,1:end,1]
            random = readdlm(plots_folder*"/"*plot_filename[1:4]*"_replicate"*string(i)*"_processed_random_net_downsampled.csv", ',')[1:end,1:end,1]
            data_centroids = replace_nan(mapslices(centroid, data, dims=1)[1,:])
            in_out_centroids = replace_nan(mapslices(centroid, in_out, dims=1)[1,:])
            out_in_centroids = replace_nan(mapslices(centroid, out_in, dims=1)[1,:])
            random_centroids = replace_nan(mapslices(centroid, random, dims=1)[1,:])
            push!(data_all, data_centroids)
            push!(in_out_all, in_out_centroids)
            push!(out_in_all, out_in_centroids)
            push!(random_all, random_centroids)
        end
        data_mean = [mean([x[i] for x in data_all]) for i in 1:length(data_all[1])]
        data_std = [std([x[i] for x in data_all]) for i in 1:length(data_all[1])]
        in_out_mean = [mean([x[i] for x in in_out_all]) for i in 1:length(in_out_all[1])]
        in_out_std = [std([x[i] for x in in_out_all]) for i in 1:length(in_out_all[1])]
        out_in_mean = [mean([x[i] for x in out_in_all]) for i in 1:length(out_in_all[1])]
        out_in_std = [std([x[i] for x in out_in_all]) for i in 1:length(out_in_all[1])]
        random_mean = [mean([x[i] for x in random_all]) for i in 1:length(random_all[1])]
        random_std = [std([x[i] for x in random_all]) for i in 1:length(random_all[1])]
        n = 5 
        ytick_interval = n/0.065/30
        data = readdlm(plots_folder*"/"*plot_filename[1:4]*"_replicate5_processed_data_intensity.csv", ',')[1:end,1:end,1]
        xs = 0:6:length(data_mean)-1
        ys = 0:ytick_interval:size(data, 1)-1
        plot_ylabel = "Distance from center (µm)"
        plot_xlabel = "Time (h)"
        fig = Figure(size=(5.5*72, 3*72))
        ax = Axis(fig[1, 1])
        colormap = Makie.wong_colors()
        lines!(ax, 1:length(data_mean), data_mean, color=:black, label="Data", linewidth=2)
        lines!(ax, 1:length(data_mean), random_mean, color=:black, linestyle=:dash, label="Random", linewidth=2)
        writedlm("$(plots_folder)/$(plot_filename[1:4])_centroids_mean.csv", data_mean, ",")
        writedlm("$(plots_folder)/$(plot_filename[1:4])_random_centroids_mean.csv", random_mean, ",")
        ax.xticks = xs
        ax.yticks = ys
        ax.xtickformat=values->string.([Int(div(v,6)) for v in values])
        ax.ytickformat=values->string.([round(Int,v*n/ytick_interval) for v in values])
        ax.xlabel = plot_xlabel
        ax.ylabel = plot_ylabel
        ax.rightspinevisible = false
        ax.topspinevisible = false
        ax.xgridvisible = false
        ax.ygridvisible = false
        fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)
        ylims!(ax, 0, nothing)
        save("$plots_folder/$plot_filename"*".pdf", fig)
    end
end

main()
