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
    images_folder = "/mnt/h/Dispersal/WT_replicate3_processed/"
    plots_folder = "/mnt/h/Dispersal/Plots"
    mask_files = sort([f for f in readdir(images_folder, join=true) if occursin("stack", f)], lt=natural)
    files = sort([f for f in readdir(folder, join=true) if occursin("piv", f)], lt=natural)

    image = load(mask_files[1])

    x, y, z, u_dummy = load(files[1], "x", "y", "z", "u")
    grid_interval =1 
    x = Int.(x)[1:grid_interval:end]
    y = Int.(y)[1:grid_interval:end]
    z = Int.(z)[1:grid_interval:end]

    u_growth = zeros(Float32, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])
    v_growth = zeros(Float32, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])
    w_growth = zeros(Float32, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])
    flags_tot = zeros(Bool, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])

    for i in 1:25
        u, v, w, flags = load(files[i], "u", "v", "w", "flags")
        u[flags .> 0] .= 0
        v[flags .> 0] .= 0
        w[flags .> 0] .= 0
        u_growth = u_growth + permutedims(u, (2,1,3))
        v_growth = v_growth + permutedims(v, (2,1,3)) 
        w_growth = w_growth + permutedims(w, (2,1,3))
        #if i < 9
        #    flags_tot += permutedims(flags, (2,1,3))
        #end
    end
    u_growth[flags_tot .> 0] .= NaN 
    v_growth[flags_tot .> 0] .= NaN  
    w_growth[flags_tot .> 0] .= NaN
    u_growth = u_growth[1:grid_interval:end, 1:grid_interval:end, 1:grid_interval:end]
    v_growth = v_growth[1:grid_interval:end, 1:grid_interval:end, 1:grid_interval:end]
    w_growth = w_growth[1:grid_interval:end, 1:grid_interval:end, 1:grid_interval:end]

    u_growth = mapwindow(median, u_growth, (3,3,1))
    v_growth = mapwindow(median, v_growth, (3,3,1))
    w_growth = mapwindow(median, w_growth, (3,3,1))
    
    distance = sqrt.(u_growth.^2 .+ v_growth.^2 .+ w_growth.^2)

    thresh=15
    u_growth[distance .< thresh] .= NaN
    v_growth[distance .< thresh] .= NaN
    w_growth[distance .< thresh] .= NaN
    
    #fig = Figure(size=(400*0.01,size(image,1)/size(image,2)*400*0.01))
    fig = Figure()
    ax = CairoMakie.Axis(fig[1,1])
    image_slice = 15 
    slice = round(Int, image_slice/size(image,3)*size(u_growth, 3))
    hm = heatmap!(ax, Int.(1:size(image, 2)).-1, Int.(size(image,1):-1:1).-1, Float32.(transpose(image[:,:,image_slice])), colormap=:grays)
    ar = arrows!(ax, (x.-1), (y.-1), u_growth[:,:,slice], v_growth[:,:,slice], arrowsize=6, arrowcolor=:salmon, linecolor=:salmon)
    writedlm("$(plots_folder)/channels_u.csv", u_growth[:,:,slice], ",")
    writedlm("$(plots_folder)/channels_v.csv", v_growth[:,:,slice], ",")
    ax.title = ""
    ax.xticksvisible=false
    ax.yticksvisible=false
    ax.xticklabelsvisible=false
    ax.yticklabelsvisible=false
    ax.rightspinevisible = false
    ax.leftspinevisible = false
    ax.topspinevisible = false
    ax.bottomspinevisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    hlines!(ax, round(Int, size(image[:,:,image_slice],2)*1/20),xmin=(round(Int, size(image[:,:,image_slice],1)*39/40) - 10/0.065)/size(image[:,:,image_slice],1),  xmax=round(Int, size(image[:,:,image_slice],1)*39/40)/size(image[:,:,image_slice],1), linewidth=10, color = :white)
    linecenter = round(Int, size(image[:,:,image_slice],1))*37/40  
    text!(ax, linecenter, round(Int, size(image[:,:,image_slice],2)*1/20)+50, text = "10 µm", align = (:left, :center), fontsize=35, color=:white)
    save(plots_folder*"/image_channel.png", fig)
end
main()
