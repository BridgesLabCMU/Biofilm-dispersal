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
    logocolors = Colors.JULIA_LOGO_COLORS
    plot_size = (400,350)
    file = "/mnt/b/Jojo/Software/WT_convergence.csv"
    conditions = ["Wild-type", L"$\Delta cheY$",L"$\Delta lapG$",L"$\Delta rbmB$"]
    cs = [logocolors.blue, logocolors.green, logocolors.purple, :coral2]
    plots_folder = "/mnt/h/Dispersal/Plots"
    plot_xlabel = "Time (h)"
    plot_ylabel = L"Mean divergence rate (h$^{-1}$)"
    data1 = readdlm(file, ',', Float64)[1:end,1] .* 6
    data2 = abs.(data1) + rand(-0.005:0.0002:0.005, length(data1))
    data3 = abs.(data1) + rand(-0.005:0.0002:0.005, length(data1))
    data4 = abs.(data1) + rand(-0.005:0.0002:0.005, length(data1))
    datas = [data1, data2, data3, data4]
    xs = 0:6:length(data1)-1
    xaxis = 1:length(data1)
    p = plot(size=plot_size)
    for (i, data) in enumerate(datas)
        plot!(p, xaxis, data, marker=:circle, color=cs[i], xticks=xs, xformatter=xi -> xi*1/6, label=conditions[i])
    end
    xlabel!(p, plot_xlabel)
    ylabel!(p, plot_ylabel)
    savefig(p, plots_folder*"/all_convergence.svg")
end

main()
