using Plots
using LaTeXStrings
using NaturalSort: sort, natural
using StatsBase

pgfplotsx()

function main()
    # Data are stored inside a folder and filename is "data.jld2" -- data are called "data"
	push!(PGFPlotsX.CUSTOM_PREAMBLE, "\\usepackage{sfmath}\n\\renewcommand{\\familydefault}{\\sfdefault}")
    default(titlefont = 20, legendfontsize = 15, 
            guidefont = (20, :black), colorbar_tickfontsize=15, colorbar_titlefontsize=20, tickfont = (15, :black), 
            guide = L"x", linewidth=2, grid=false, formatter=:plain)
    plot_xlabel = "Time (h)"
    plot_ylabel = "Reporter activity (a.u.)"

    master_directory = "/mnt/f/Sandhya_Imaging/Time_Lapses/Reporters/"
    data_folders = sort(filter(x -> occursin("data", x), readdir(master_directory, join=true)), lt=natural) 
    plots_folder = "/mnt/h/Dispersal/Plots"
    plot_size = (350,300*length(plot_folders))
    p = plot(layout=(1,length(data_folders)), size=plot_size)
    titles = [L"\textit{Ptac}", L"\textit{vpsL}", L"\textit{flaA}", L"\textit{luxC}"]

    for i in length(data_folders)
        data = load(data_folder*"data.jld2", "data")
        x=1:length(data)
        plot!(x, data, xlabel=plot_xlabel, ylabel=plot_ylabel, title=titles[i], xformatter=xi -> xi*1/2, subplot=i)
    end
	savefig("$plots_folder/$plot_filename"*".svg")
end
