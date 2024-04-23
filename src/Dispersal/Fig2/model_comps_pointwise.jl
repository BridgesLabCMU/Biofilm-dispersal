using Plots, StatsPlots
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
using Distances

pgfplotsx()

function main()
    push!(PGFPlotsX.CUSTOM_PREAMBLE, "\\usepackage{sfmath}\n\\renewcommand{\\familydefault}{\\sfdefault}")
    default(titlefont = 20, legendfontsize = 15, 
            guidefont = (20, :black), colorbar_tickfontsize=15, colorbar_titlefontsize=20, tickfont = (15, :black), 
            guide = L"x", linewidth=2, grid=false, formatter=:plain)
    plot_size = (300,250)
    plot_size2 = (300,250)
    plot_xlabel = "Time (h)"
    plot_ylabel = "J-S divergence"
    plots_folder = "/mnt/h/Dispersal/Plots"
    files = [f for f in readdir(plots_folder, join=true) if occursin("data", f) && occursin("csv", f) && occursin("WT", f)]
    in_out_total = Array{Float64, 1}(undef, length(files))
    out_in_total = Array{Float64, 1}(undef, length(files))
    random_total = Array{Float64, 1}(undef, length(files))
    colors = theme_palette(:auto)
    for (i, file) in enumerate(files)
        data = readdlm(file, ',')[1:end,1:end,1]
        in_out = readdlm(file[1:end-9]*"_in_out_pointwise.csv", ',')[1:end,1:end,1]
        out_in = readdlm(file[1:end-9]*"_out_in_pointwise.csv", ',')[1:end,1:end,1]
        random = readdlm(file[1:end-9]*"_random_pointwise.csv", ',')[1:end,1:end,1]
        in_out_comp = Array{Float64}(undef, size(data,2))
        out_in_comp = Array{Float64}(undef, size(data,2))
        random_comp = Array{Float64}(undef, size(data,2))
        for j in 1:size(data, 2)
            data_norm = data[:,j] ./ sum(data[:,j])
            in_out_comp[j] = js_divergence(in_out[:,j] ./ sum(in_out[:,j]), data_norm)
            out_in_comp[j] = js_divergence(out_in[:,j] ./ sum(out_in[:,j]), data_norm)
            random_comp[j] = js_divergence(random[:,j] ./ sum(random[:,j]), data_norm)
        end
        in_out_total[i] = mean(filter(!isnan, in_out_comp)) 
        out_in_total[i] = mean(filter(!isnan, out_in_comp))
        random_total[i] = mean(filter(!isnan, random_comp))
        xs = 0:6:size(data,2)-1
        xaxis = 1:size(data,2)
        p = plot(ylim=(0,1))
        scatter!(p, xaxis, in_out_comp, xticks=xs, xformatter=xi -> xi*1/6, size=plot_size, label="Inside-out")
        scatter!(p, xaxis, out_in_comp, xticks=xs, xformatter=xi -> xi*1/6, size=plot_size, label="Outside-in")
        scatter!(p, xaxis, random_comp, xticks=xs, xformatter=xi -> xi*1/6, size=plot_size, label="Random")
        xlabel!(p, plot_xlabel)
        ylabel!(p, plot_ylabel)
        savefig(p, file[1:end-9]*"_comps_pointwise.svg")
    end
    p = plot(ylim=(0,1), xticks=([1,2,3],["Inside-out", "Outside-in", "Random"]), xrotation=45, size=plot_size2, leg=false)
    boxplot!(p, [in_out_total, out_in_total, random_total], linewidth = 1.0, outliers=false)
    xlabel!(p, "")
    ylabel!(p, "Average J-S divergence")
    savefig(p, files[1][1:end-18]*"_comps_pointwise_summarized.svg")
end

main()
