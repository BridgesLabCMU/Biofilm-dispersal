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
using Interpolations
using FileIO

pgfplotsx()

round_up(x, multiple) = ceil(x / multiple) * multiple

function deform_image!(deformed_image, image, displacements_folder, first_timepoint, current_timepoint)
    x_grid, y_grid, z_grid, u_tot = load("$displacements_folder/piv_results_$(first_timepoint).jld2", 
                                         "x", "y", "z", "u")
    y_grid = y_grid[end:-1:1]
    u_tot .= 0
    v_tot = zeros(Float32, size(u))
    w_tot = zeros(Float32, size(u))
    for t in first_timepoint:current_timepoint-1
        u, v, w = load("$displacements_folder/piv_results_$(t).jld2", "u", "v", "w")
        u_tot .+= u
        v_tot .+= v
        w_tot .+= w
    end
    v_tot .*= -1
    h, w, d = size(image)
    x = 1:w
    y = 1:h
    z = 1:d
    itp_u_img = extrapolate(scale(interpolate(u_tot, BSpline(Cubic(Line(OnGrid())))), 
                                  (y_grid, x_grid, z_grid)), Line())
    itp_v_img = extrapolate(scale(interpolate(v_tot, BSpline(Cubic(Line(OnGrid())))), 
                                  (y_grid, x_grid, z_grid)), Line())
    itp_w_img = extrapolate(scale(interpolate(w_tot, BSpline(Cubic(Line(OnGrid())))), 
                                  (y_grid, x_grid, z_grid)), Line())
    ut = itp_u_img(y, x, z)
    vt = itp_v_img(y, x, z)
    wt = itp_w_img(y, x, z)
    itp_img = extrapolate(interpolate(image, BSpline(Cubic(Line(OnGrid())))), Line())
    deformed_image .= itp_img.(y .+ vt, x .+ ut, z .+ wt)
    return nothing
end

function mask_image!(intensity_thresholds, mask, deformed_image)
    slices = size(deformed_image, 3)
    for i in 1:slices
        masks[:,:,i] = @views deformed_image[:,:,i] .> intensity_thresholds[i]
    end
end

function radial_averaging(mask_files, files, images_folder, first_index, end_index, bin_interval, intensity_thresholds, displacements_folder)
    first_image = load("$images_folder/$(mask_files[first_index])")
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
    ntimepoints = end_index - first_index + 1
    data_matrix = zeros(nbins-1, ntimepoints)
    for t in 1:ntimepoints
        image = load("$images_folder/$(files[t+first_index-1])")
        mask = zeros(Bool, size(image))
        if t == 1
            mask_image!(intensity_thresholds, mask, image)
        else
            deformed_image = similar(image)
            deform_image!(deformed_image, image, displacements_folder, first_index, t+first_index-1)
            mask_image!(intensity_thresholds, mask, deformed_image)
        end
        for i in 1:size(data_matrix, 1) 
            data_matrix[i, t] = mean(mask[findall(x -> bins[i] <= x <= bins[i+1], center_distance)])
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
    plot_title = "Data"

    master_directory = "/mnt/h/Dispersal"
    image_folders = filter(isdir, readdir(master_directory, join=true))
    filter!(folder->folder≠master_directory*"/Plots", image_folders)
    plots_folder = "/mnt/h/Dispersal/Plots"

    for images_folder in image_folders
        plot_filename = basename(images_folder)*"_data" 
        displacements_folder = "$(images_folder)/Displacements"
        intensity_thresholds_file = "$(images_folder)/intensity_thresholds.csv" 
        files = sort([f for f in readdir(images_folder) if occursin("noplank_isotropic", f)], 
                                 lt=natural)
        files = sort([f for f in readdir(images_folder) if occursin("mask_isotropic", f)], 
                                 lt=natural)
        intensity_thresholds = readdlm(intensity_thresholds_folder, ',', Float64)[1:end,1]
        ntimepoints = length(files)
        first_image = load("$images_folder/$(files[1])"; lazyio=true)
        height, width, slices = size(first_image)
        first_image = nothing
        net = readdlm(plots_folder*"/"*basename(images_folder)*".csv", ',', Int)[1:end,1]
        first_index = argmax(net)
        end_index = min(first_index + 45, ntimepoints)

        data_matrix = diff(radial_averaging(mask_files, files, images_folder, first_index, end_index, 30, intensity_thresholds, displacements_folder), dims=2)
        data_matrix[data_matrix .> 0] .= 0
        replace!(data_matrix, Inf=>NaN)
        replace!(data_matrix, NaN=>0.0)
        writedlm("$(plots_folder)/$(plot_filename).csv", data_matrix, ",")
        n = 10
        ytick_interval = n/0.065/30
        xs = 0:6:size(data_matrix, 2)-1
        ys = 0:ytick_interval:size(data_matrix, 1)-1
        c = cgrad(:Purples, rev=true)
        plt = heatmap(data_matrix, xticks=xs, yticks=ys, color=c, clim=(-0.1, 0),
                      colorbar_title="Density change (a.u.)", xformatter=xi -> xi*1/6, 
                      yformatter=yi -> yi/ytick_interval*n, size=plot_size)
        xlabel!(plt, plot_xlabel)
        ylabel!(plt, plot_ylabel)
        title!(plt, plot_title)
        savefig("$plots_folder/$plot_filename"*".svg")
    end
end

main()
