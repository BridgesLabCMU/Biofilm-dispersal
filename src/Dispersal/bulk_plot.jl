using Plots
using LaTeXStrings
using StatsBase
using DelimitedFiles
using Colors

pgfplotsx()

function main()
    push!(PGFPlotsX.CUSTOM_PREAMBLE, "\\usepackage{sfmath}\n\\renewcommand{\\familydefault}{\\sfdefault}")
    default(titlefont = 20, legendfontsize = 15, 
            guidefont = (20, :black), colorbar_tickfontsize=15, colorbar_titlefontsize=20, tickfont = (15, :black), 
            guide = L"x", linewidth=2, grid=false, formatter=:plain)
    plot_size = (400,350)

    plots_folder = "/mnt/h/Dispersal/Plots"
    files = [f for f in readdir(plots_folder, join=true) if occursin(".csv", f)]
    plot_xlabel = "Time (h)"
    plot_ylabel = "Biovolume (a.u.)" 
    logocolors = Colors.JULIA_LOGO_COLORS
    p = plot(ylim=(0,1))
    WT_seen = false
    cheY_seen = false
    rbmB_seen = false
    lapG_seen = false
    for file in files
        data = readdlm(file, ',', Int)[1:end,1]
        first_idx = argmax(data)
        end_idx = min(first_idx + 45, length(data))
        data_subset = data[first_idx:end_idx] 
        data = data_subset
        data_norm =  data ./ data[1]
        data = data_norm
        lab = split(basename(file), "_")[1]
        if lab == "WT"
            c = logocolors.blue
            if WT_seen
                condition = ""
            else
                condition = "Wild-type"
                WT_seen = true
            end
        elseif lab == "cheY"
            c = logocolors.green
            if cheY_seen
                condition = ""
            else
                condition = latexstring("\\Delta $(lab)")
                cheY_seen = true
            end
        elseif lab == "rbmB"
            c = :coral2 
            if rbmB_seen
                condition = ""
            else
                condition = latexstring("\\Delta $(lab)")
                rbmB_seen = true
            end
        elseif lab == "lapG"
            c = logocolors.purple
            if lapG_seen
                condition = ""
            else
                condition = latexstring("\\Delta $(lab)")
                lapG_seen = true
            end
        end
        xs = 0:6:length(data)-1
        xaxis = 1:length(data)
        plot!(p, xaxis, data, marker=:circle, color=c, xticks=xs, xformatter=xi -> xi*1/6, label=condition)
    end
    xlabel!(p, plot_xlabel)
    ylabel!(p, plot_ylabel)
    savefig(p, plots_folder*"/WT_cheY_lapG_rbmB_bulk.svg")
end

main()
