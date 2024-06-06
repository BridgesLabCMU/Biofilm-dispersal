using FileIO 
using Makie 
using GLMakie
using CairoMakie
CairoMakie.activate!(type = "pdf")
set_theme!(fonts = (; regular = "TeX Gyre Heros Regular", bold = "TeX Gyre Heros Regular"))
using NaturalSort
using StatsBase
using NaNStatistics
using ImageMorphology
using Images

function main()
    folder = "/mnt/h/Dispersal/WT_replicate3_processed/Displacements/"
    plots_folder = "/mnt/h/Dispersal/Plots"
    files = sort([f for f in readdir(folder, join=true) if occursin("piv", f)], lt=natural)

    x, y, z, u_dummy = load(files[1], "x", "y", "z", "u")
    x = Int.(x)
    y = Int.(y)
    z = Int.(z)

    flags_tot = zeros(Bool, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])

    n_out = Array{Float32, 1}(undef, length(files))

    for i in 1:length(files)
        u, v, w, flags = load(files[i], "u", "v", "w", "flags")
        u[flags .> 0] .= 0
        v[flags .> 0] .= 0
        w[flags .> 0] .= 0
        #u_growth = u_growth + permutedims(u, (2,1,3))
        #v_growth = v_growth + permutedims(v, (2,1,3)) 
        #w_growth = w_growth + permutedims(w, (2,1,3))
        u_growth = mapwindow(median, u, (3,3,1))
        v_growth = mapwindow(median, v, (3,3,1))
        w_growth = mapwindow(median, w, (3,3,1))
        u_mean = mapwindow(mean, u_growth, (3,3,1)) 
        v_mean = mapwindow(mean, v_growth, (3,3,1)) 
        w_mean = mapwindow(mean, w_growth, (3,3,1)) 
        u_out = u_growth .- u_mean
        v_out = v_growth .- v_mean
        w_out = w_growth .- w_mean
        distance = sqrt.(u_out.^2 .+ v_growth.^2 .+ w_growth.^2)
        thresh=15
        n_out[i] = sum(distance .> thresh)
        #u_growth[distance .> thresh] .= NaN
        #v_growth[distance .> thresh] .= NaN
        #w_growth[distance .> thresh] .= NaN
    end
    
    n_out .*= 8^3
    n_out .*= 0.065^3
    fig = Figure()
    ax = CairoMakie.Axis(fig[1,1])
    ax.title = ""
    xs = 0:6:length(files)-1
    xaxis = 1:length(files)
    xaxis /= 6
    lines!(ax, xaxis, n_out, color=:black, alpha=0.3)
    ax.xlabel = "Time (h)"
    ax.ylabel = rich("Channel volume (μm", superscript("3"), ")")
    ax.rightspinevisible = false
    ax.topspinevisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    save(plots_folder*"/WT_channels.pdf", fig)
end
main()
