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
	plt = heatmap(data_matrix, xticks=xs, yticks=ys, color=c, clim=(-0.1, 0),
				  colorbar_title="Density change (a.u.)", xformatter=xi -> xi*1/6, 
				  yformatter=yi -> yi/ytick_interval*n, size=plot_size)
	xlabel!(plt, plot_xlabel)
	ylabel!(plt, plot_ylabel)
	title!(plt, plot_title)
	savefig("$plots_folder/$plot_filename"*".svg")
end
