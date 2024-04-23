using Plots
using LaTeXStrings
using StatsBase
using DelimitedFiles
using Colors
using NaturalSort
using FileIO

pgfplotsx()

function choose_label(folder)
    if occursin("WT", folder)
        return "WT"
    elseif occursin("cheY", folder)
        return "cheY"
    elseif occursin("lapG", folder)
        return "lapG"
    elseif occursin("rbmB", folder)
        return "rbmB"
    end
end

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
    net_divergence /= nx*ny*nz
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
    p = plot(size=plot_size)
    vector_folders = ["/mnt/h/Dispersal/WT_replicate1_processed/Displacements/", 
                      "/mnt/h/Dispersal/WT_replicate2_processed/Displacements/", 
                      "/mnt/h/Dispersal/WT_replicate3_processed/Displacements/", 
                      "/mnt/h/Dispersal/WT_replicate4_processed/Displacements/", 
                      "/mnt/h/Dispersal/WT_replicate5_processed/Displacements/", 
                      "/mnt/h/Dispersal/cheY_replicate1_processed/Displacements/", 
                      "/mnt/h/Dispersal/cheY_replicate2_processed/Displacements/", 
                      "/mnt/h/Dispersal/cheY_replicate3_processed/Displacements/", 
                      "/mnt/h/Dispersal/cheY_replicate4_processed/Displacements/", 
                      "/mnt/h/Dispersal/cheY_replicate5_processed/Displacements/", 
                      "/mnt/h/Dispersal/lapG_replicate1_processed/Displacements/", 
                      "/mnt/h/Dispersal/lapG_replicate2_processed/Displacements/", 
                      "/mnt/h/Dispersal/lapG_replicate3_processed/Displacements/", 
                      "/mnt/h/Dispersal/lapG_replicate4_processed/Displacements/", 
                      "/mnt/h/Dispersal/lapG_replicate5_processed/Displacements/", 
                      "/mnt/h/Dispersal/rbmB_replicate1_processed/Displacements/", 
                      "/mnt/h/Dispersal/rbmB_replicate2_processed/Displacements/", 
                      "/mnt/h/Dispersal/rbmB_replicate3_processed/Displacements/", 
                      "/mnt/h/Dispersal/rbmB_replicate4_processed/Displacements/", 
                      "/mnt/h/Dispersal/rbmB_replicate5_processed/Displacements/"]
    dx = 8
    dy = 8
    dz = 8
    WT_seen = false
    cheY_seen = false
    rbmB_seen = false
    lapG_seen = false
    bulk_files = [f for f in readdir(plots_folder, join=true) if occursin("processed.csv", f)]
    logocolors = Colors.JULIA_LOGO_COLORS
    WT_averages = []
    cheY_averages = []
    rbmB_averages = []
    lapG_averages = []
    for vector_folder in vector_folders
        bulk_file = [f for f in bulk_files if occursin(split(vector_folder, "/")[end-2], f)][1]
        bulk_data = readdlm(bulk_file, ',', Int)[1:end,1]
        first_idx = argmax(bulk_data)
        end_idx = min(first_idx+45, length(bulk_data))
        rpca_files = [f for f in readdir(vector_folder, join=true) if occursin("rpca", f)]
        u, v, w = load(rpca_files[1], "u", "v", "w")
        divergences = Array{Float64, 1}(undef, end_idx-first_idx+1)
        for i in first_idx:end_idx 
            @views ui = permutedims(u[:,:,:,i], [2,1,3])
            @views vi = permutedims(v[:,:,:,i], [2,1,3])
            @views wi = permutedims(w[:,:,:,i], [2,1,3])
            net_divergence = calculate_divergence(ui, vi, wi, dx, dy, dz)
            divergences[i-first_idx+1] = net_divergence
        end
        xs = 0:6:length(divergences)-1
        xaxis = 1:length(divergences)
        lab = choose_label(vector_folder)
        if lab == "WT"
            push!(WT_averages, divergences)
            c = logocolors.blue
            if WT_seen
                condition = ""
            else
                condition = "Wild-type"
                WT_seen = true
            end
        elseif lab == "cheY"
            push!(cheY_averages, divergences)
            c = logocolors.green
            if cheY_seen
                condition = ""
            else
                condition = latexstring("\\Delta $(lab)")
                cheY_seen = true
            end
        elseif lab == "rbmB"
            push!(rbmB_averages, divergences)
            c = :coral2 
            if rbmB_seen
                condition = ""
            else
                condition = latexstring("\\Delta $(lab)")
                rbmB_seen = true
            end
        elseif lab == "lapG"
            push!(lapG_averages, divergences)
            c = logocolors.purple
            if lapG_seen
                condition = ""
            else
                condition = latexstring("\\Delta $(lab)")
                lapG_seen = true
            end
        end
        scatter!(p, xaxis, divergences, mc=c, ma=0.1, xticks=xs, xformatter=xi -> xi*1/6, label="")
    end
    WT_averages = hcat(WT_averages...)
    WT_averages = mean(WT_averages, dims=2)
    cheY_averages = hcat(cheY_averages...)
    cheY_averages = mean(cheY_averages, dims=2)
    rbmB_averages = hcat(rbmB_averages...)
    rbmB_averages = mean(rbmB_averages, dims=2)
    lapG_averages = hcat(lapG_averages...)
    lapG_averages = mean(lapG_averages, dims=2)
    xaxis = 1:length(WT_averages)
    xs = 0:6:length(WT_averages)-1
    plot!(p, xaxis, WT_averages, c=logocolors.blue, xticks=xs, xformatter=xi -> xi*1/6, label="Wild-type")
    plot!(p, xaxis, cheY_averages, c=logocolors.green, xticks=xs, xformatter=xi -> xi*1/6, label=L"$\Delta cheY$")
    plot!(p, xaxis, rbmB_averages, c=:coral2, xticks=xs, xformatter=xi -> xi*1/6, label=L"$\Delta rbmB$")
    plot!(p, xaxis, lapG_averages, c=logocolors.purple, xticks=xs, xformatter=xi -> xi*1/6, label=L"$\Delta lapG$")
    xlabel!(p, plot_xlabel)
    ylabel!(p, plot_ylabel)
    savefig(p, plots_folder*"/all_convergence.svg")
end

main()
