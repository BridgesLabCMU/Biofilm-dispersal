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

pgfplotsx()

function model_comparison(A, B)
    # A = Data, B = model
    net_A = sum(A)
    net_B = sum(B)
    D = 0
    for i in eachindex(A)
        D += A[i] * (log(A[i]/B[i]) + log(net_B) - log(net_A))
    end
    return D
end

function main()
    push!(PGFPlotsX.CUSTOM_PREAMBLE, "\\usepackage{sfmath}\n\\renewcommand{\\familydefault}{\\sfdefault}")
    default(titlefont = 20, legendfontsize = 15, 
            guidefont = (20, :black), colorbar_tickfontsize=15, colorbar_titlefontsize=20, tickfont = (15, :black), 
            guide = L"x", linewidth=2, grid=false, formatter=:plain)
    plot_size2 = (300,250)
    plots_folder = "/mnt/h/Dispersal/Plots"
    files = [f for f in readdir(plots_folder, join=true) if occursin("data", f) && occursin("csv", f) && occursin("lapG", f)]
    @show files
    in_out_total = Array{Float64, 1}(undef, length(files))
    out_in_total = Array{Float64, 1}(undef, length(files))
    random_total = Array{Float64, 1}(undef, length(files))
    colors = theme_palette(:auto)
    for (i, file) in enumerate(files)
        data = readdlm(file, ',')[1:end,1:end,1]
        in_out = readdlm(file[1:end-9]*"_in_out.csv", ',')[1:end,1:end,1]
        out_in = readdlm(file[1:end-9]*"_out_in.csv", ',')[1:end,1:end,1]
        random = readdlm(file[1:end-9]*"_random.csv", ',')[1:end,1:end,1]
        in_out_total[i] = model_comparison(data, in_out)
        out_in_total[i] = model_comparison(data, out_in)
        random_total[i] = model_comparison(data, random)
    end
    p = plot(ylim=(0,1), xticks=([1,2,3],["Inside-out", "Outside-in", "Random"]), xrotation=45, size=plot_size2, leg=false)
    boxplot!(p, [in_out_total, out_in_total, random_total], linewidth = 1.0, outliers=false)
    xlabel!(p, "")
    ylabel!(p, "Similarity (a.u.)")
    savefig(p, files[1][1:end-18]*"_comps_summarized.svg")
end

main()
