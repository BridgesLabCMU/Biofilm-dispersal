using FileIO 
using NPZ
using CairoMakie
using NaturalSort
using StatsBase
using NaNStatistics

#x, y, z, u, v, w, s2n, flags = load("../tests/Displacements/piv_results_2.jld2", 
#                              "x", "y", "z", "u", "v", "w", "s2n", "flags")
#x_p = npzread("../tests/Displacements/x.npy")
#y_p = npzread("../tests/Displacements/y.npy")
#z_p = npzread("../tests/Displacements/z.npy")
#u_p = npzread("../tests/Displacements/u.npy")
#v_p = npzread("../tests/Displacements/v.npy")
#w_p = npzread("../tests/Displacements/w.npy")
#s2n_p = npzread("../tests/Displacements/s2n.npy")
#flags_p = npzread("../tests/Displacements/flags.npy") 

#flags_qplot = transpose(flags[:,:,8])
#u_qplot = transpose(u[:,:,8])
#v_qplot = transpose(v[:,:,8])
#u_qplot[flags_qplot] .= 0.0
#v_qplot[flags_qplot] .= 0.0

#u_p = transpose(u_p[2,:,:])
#v_p = transpose(v_p[2,:,:])
#flags_p = transpose(flags_p[2,:,:])
#u_p[flags_p] .= 0.0  
#v_p[flags_p] .= 0.0  

# 2D quiver plot of u_qplot, v_qplot
#fig = arrows((x .- 1) ./ 8, (y .- 1) ./ 8, u_p, v_p)
#save("../tests/Displacements/quiver_plot_python_postinterp.png", fig)

function main()
    folder = "/mnt/h/Dispersal/WT_replicate1_processed/Displacements/"
    rpca_file = [f for f in readdir(folder, join=true) if occursin("rpca", f)]
    files = 

    x, y, z, u_dummy = load(files[1], "x", "y", "z", "u")
    x = Int.(x)
    y = Int.(y)
    z = Int.(z)

    u_p = zeros(Float32, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])
    v_p = zeros(Float32, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])
    w_p = zeros(Float32, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])
    flags_tot = zeros(Bool, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])

    for i in 1:10
        u, v, w, flags = load(files[i], "u", "v", "w", "flags")
        u_p += permutedims(u, (2,1,3))
        v_p += permutedims(v, (2,1,3)) 
        w_p += permutedims(w, (2,1,3))
        flags_tot += permutedims(flags, (2,1,3))
    end

    u_p[flags_tot .> 0] .= NaN 
    v_p[flags_tot .> 0] .= NaN  
    w_p[flags_tot .> 0] .= NaN

    #u_p[distance .< 5] .= NaN
    #v_p[distance .< 5] .= NaN
    #w_p[distance .< 5] .= NaN

    #ps = [Point3f(xi/8, yi/8, zi/2) for xi in x for yi in y for zi in z]
    #ns = [Vec3f(u_p[i], v_p[i], w_p[i]) for i in 1:length(u_p)]
    fig = arrows((x.-1), (z.-1).*4, u_p[:,30,:], w_p[:,30,:].*4)
    #fig = arrows(ps, ns, fxaa=true)
    save("$folder/quiver_rpca.png", fig)
end
main()
