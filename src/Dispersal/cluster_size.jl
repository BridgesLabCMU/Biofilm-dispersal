using Plots
using TiffImages: load, save
using ImageMorphology: label_components, component_lengths
using StatsBase: Histogram, fit
using Interpolations
using FileIO

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

function cluster_sizes!(binary_timeseries, ntimepoints, conversion_factor, clusters_time)
    @inbounds for t in 2:ntimepoints
        dispersed = binary_timeseries[:,:,:,t-1] .& .~binary_timeseries[:,:,:,t]
        clusters = label_components(dispersed)
        cluster_sizes = component_lengths(clusters)[1:end] .* conversion_factor
        clusters_time[t] = cluster_sizes
    end
    return clusters_time
end

function read_images!(directory, ntimepoints, arr, files)
    @inbounds for t in 1:ntimepoints
        @views file = files[t]
        arr[:, :, :, t] = load("$directory/$file")
    end
    return nothing 
end

function sizes_to_hist(clusters_time, ntimepoints)
    flat_array = vcat(clusters_time...)
    max_size = ceil(maximum(flat_array))
    sizes = 0:1:max_size-1
    data_matrix = Array{Int, 2}(undef, length(sizes), ntimepoints)
    for t in 1:ntimepoints
        data_matrix[:, t] = fit(Histogram, clusters_time[t], sizes)
    end
    return data_matrix
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
    plot_filename = "cheY_replicate1_clusters"
    plot_xlabel = "Time (h)"
    plot_ylabel = L"Cluster size ($\mu$m$^3$)"
    conversion_factor = 0.065^3
    files = sort([f for f in readdir(images_folder) if occursin("mask", f)], 
                             lt=natural)
    ntimepoints = length(files)
    dummy_image = load("$images_folder/$(files[1])"; lazyio=true)
    height, width, slices = size(dummy_image)
    timeseries = zeros(Bool, height, width, slices, ntimepoints)
    read_images!(images_folder, ntimepoints, timeseries, files)
    clusters_time = Array{Float64, 1}(undef, ntimepoints)
    cluster_sizes!(timeseries, ntimepoints, conversion_factor, clusters_time)
    data_matrix = sizes_to_hist(clusters_time, ntimepoints)
	plt = heatmap(plot_xticks, reverse(plot_yticks), data_matrix, 
                  yflip=true, color=:coolwarm,
                  colorbar_title=clab_title, size=plot_size)
    xlabel!(plt, plot_xlabel)
    ylabel!(plt, plot_ylabel)
    title!(plt, plot_title)
    savefig("$plots_folder/$plot_filename"*".svg")
end

main()
