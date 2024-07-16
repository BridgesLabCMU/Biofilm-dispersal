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
    plot_filename = "rbmA_centroids"
    n = 5
    ytick_interval = n/0.065/30
    plot_ylabel = "Distance from center (µm)"
    plot_xlabel = "Time (h)"
    for i in 1:5
        data = readdlm(plots_folder*"/rbmA_replicate"*string(i)*"_processed_data_intensity.csv", ',')[1:end,1:end,1]
        random = readdlm(plots_folder*"/rbmA_replicate"*string(i)*"_processed_random_net_downsampled.csv", ',')[1:end,1:end,1]
        data_centroids = replace_nan(mapslices(centroid, data, dims=1)[1,:])
        random_centroids = replace_nan(mapslices(centroid, random, dims=1)[1,:])
        xs = 0:6:length(data_centroids)-1
        ys = 0:ytick_interval:size(data, 1)-1
        fig = Figure(size=(5.5*72, 3*72))
        ax = Axis(fig[1, 1])
        lines!(ax, 1:length(data_centroids), data_centroids, label="Data", color=:black, linewidth=2)
        lines!(ax, 1:length(data_centroids), random_centroids, label="Random", color=:black, linestyle=:dash, colorrange=(1,7), linewidth=2)
        writedlm("$(plots_folder)/rbmA_centroids_$(string(i)).csv", data_centroids, ",")                                                                 
        writedlm("$(plots_folder)/rbmA_random_centroids_$(string(i)).csv", random_centroids, ",")
        ax.xticks = xs
        ax.yticks = ys
        ax.xtickformat=values->string.([Int(div(v,6)) for v in values])
        ax.ytickformat=values->string.([round(Int,(v*n/ytick_interval)) for v in values])
        ax.xlabel = plot_xlabel
        ax.ylabel = plot_ylabel
        ax.rightspinevisible = false
        ax.topspinevisible = false
        ax.xgridvisible = false
        ax.ygridvisible = false
        ylims!(ax, 0, nothing)
        fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)
        save("$plots_folder/$plot_filename"*string(i)*".pdf", fig)
    end
end

main()
