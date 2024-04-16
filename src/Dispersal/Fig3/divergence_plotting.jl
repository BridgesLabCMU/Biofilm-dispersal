using Plots
using LaTeXStrings
using StatsBase
using DelimitedFiles
using Colors
using NaturalSort
using FileIO

pgfplotsx()

function calculate_divergence(u, v, w, dx, dy, dz)
    nx, ny, nz = size(u)
    net_divergence = 0.0
    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
				dudx = i==nx ? (u[i,j,k]-u[i-1,j,k])/dx : i==1 ? (u[i+1,j,k]-u[i,j,k])/dx : (u[i+1,j,k]-u[i-1,j,k])/(2dx)
                dvdy = j==ny ? (v[i,j,k]-v[i,j-1,k])/dy : j==1 ? (v[i,j+1,k]-v[i,j,k])/dy : (v[i,j+1,k]-v[i,j-1,k])/(2dy)
                dwdz = k==nz ? (w[i,j,k]-w[i,j,k-1])/dz : k==1 ? (w[i,j,k+1]-w[i,j,k])/dz : (w[i,j,k+1]-w[i,j,k-1])/(2dz)
                net_divergence += dudx + dvdy + dwdz
            end
        end
    end
    return net_divergence
end

function main()
    push!(PGFPlotsX.CUSTOM_PREAMBLE, "\\usepackage{sfmath}\n\\renewcommand{\\familydefault}{\\sfdefault}")
    default(titlefont = 20, legendfontsize = 15, 
            guidefont = (20, :black), colorbar_tickfontsize=15, colorbar_titlefontsize=20, tickfont = (15, :black), 
            guide = L"x", linewidth=2, grid=false, formatter=:plain)
    logocolors = Colors.JULIA_LOGO_COLORS
    plot_size = (400,350)
    conditions = ["Wild-type", L"$\Delta cheY$",L"$\Delta lapG$",L"$\Delta rbmB$"]
    cs = [logocolors.blue, logocolors.green, logocolors.purple, :coral2]
    plots_folder = "/mnt/h/Dispersal/Plots"
    plot_xlabel = "Time (h)"
    plot_ylabel = L"Net divergence rate (h$^{-1}$)"
    vector_folder = "/mnt/h/Dispersal/WT_replicate2_processed/Displacements/"
    dx = 8
    dy = 8
    dz = 8
    rpca_files = [f for f in readdir(vector_folder, join=true) if occursin("rpca", f)]
    u, v, w = load(rpca_files[1], "u", "v", "w")
    divergences = Array{Float64, 1}(undef, size(u)[4])
    for i in 1:size(u)[4] 
        @views ui = permutedims(u[:,:,:,i], [2,1,3])
        @views vi = permutedims(v[:,:,:,i], [2,1,3])
        @views wi = permutedims(w[:,:,:,i], [2,1,3])
        net_divergence = calculate_divergence(ui, vi, wi, dx, dy, dz)
        divergences[i] = net_divergence
    end
    xs = 0:6:length(divergences)-1
    xaxis = 1:length(divergences)
    p = plot(size=plot_size)
    plot!(p, xaxis, divergences, marker=:circle, color=logocolors.blue, xticks=xs, xformatter=xi -> xi*1/6, label="Wild-type")

    """
    datas = [data1, data2, data3, data4]
    xs = 0:6:length(data1)-1
    xaxis = 1:length(data1)
    p = plot(size=plot_size)
    for (i, data) in enumerate(datas)
        plot!(p, xaxis, data, marker=:circle, color=cs[i], xticks=xs, xformatter=xi -> xi*1/6, label=conditions[i])
    end
    """
    xlabel!(p, plot_xlabel)
    ylabel!(p, plot_ylabel)
    savefig(p, plots_folder*"/all_convergence.svg")
end

main()
