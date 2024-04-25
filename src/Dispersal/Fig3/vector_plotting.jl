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

function calculate_divergence(u, v, w, dx, dy, dz)
    nx, ny, nz = size(u)
    divergence = zeros(Float32, nx, ny, nz)
    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                if !isnan(u[i,j,k])
                    dudx = (i==nx && !isnan(u[i-1,j,k])) ? (u[i,j,k]-u[i-1,j,k])/dx : (i==1 && !isnan(u[i+1,j,k])) ? (u[i+1,j,k]-u[i,j,k])/dx : (i!=nx && i!=1 && !isnan(u[i-1,j,k]) && !isnan(u[i+1,j,k])) ? (u[i+1,j,k]-u[i-1,j,k])/(2dx) : NaN 
				else
                    dudx = NaN 
                end
                if !isnan(v[i,j,k])
                    dvdy = (j==ny && !isnan(v[i,j-1,k])) ? (v[i,j,k]-v[i,j-1,k])/dx : (j==1 && !isnan(v[i,j+1,k])) ? (v[i,j+1,k]-v[i,j,k])/dx : (j!=ny && j!=1 && !isnan(v[i,j-1,k]) && !isnan(v[i,j+1,k])) ? (v[i,j+1,k]-v[i,j-1,k])/(2dx) : NaN 
                else
                    dvdy = NaN 
                end
                if !isnan(w[i,j,k])
                    dwdz = (k==nz && !isnan(w[i,j,k-1])) ? (w[i,j,k]-w[i,j,k-1])/dx : (k==1 && !isnan(w[i,j,k+1])) ? (w[i,j,k+1]-w[i,j,k])/dx : (k!=nz && k!=1 && !isnan(w[i,j,k-1]) && !isnan(w[i,j,k+1])) ? (w[i,j,k+1]-w[i,j,k-1])/(2dx) : NaN 
                else
                    dwdz = NaN 
                end
                divergence[i,j,k] = dudx + dvdy + dwdz
            end
        end
    end
    return divergence
end

function main()
    folder = "/mnt/h/Dispersal/WT_replicate1_processed/Displacements/"
    #rpca_file = [f for f in readdir(folder, join=true) if occursin("rpca", f)]

    files = sort([f for f in readdir(folder, join=true) if occursin("piv", f)], lt=natural)

    x, y, z, u_dummy = load(files[1], "x", "y", "z", "u")
    x = Int.(x)
    y = Int.(y)
    z = Int.(z)

    u_p = zeros(Float32, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])
    v_p = zeros(Float32, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])
    w_p = zeros(Float32, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])
    flags_tot = zeros(Bool, size(u_dummy)[2], size(u_dummy)[1], size(u_dummy)[3])

    for i in 1:15 
        u, v, w, flags = load(files[i], "u", "v", "w", "flags")
        u_p += permutedims(u, (2,1,3))
        v_p += permutedims(v, (2,1,3)) 
        w_p += permutedims(w, (2,1,3))
        flags_tot += permutedims(flags, (2,1,3))
    end

    #u_p[flags_tot .> 0] .= NaN 
    #v_p[flags_tot .> 0] .= NaN  
    #w_p[flags_tot .> 0] .= NaN
    u_p[flags_tot .> 0] .= NaN 
    v_p[flags_tot .> 0] .= NaN  
    w_p[flags_tot .> 0] .= NaN
    divergence = calculate_divergence(u_p, v_p, w_p.*4, 8, 8, 8)

    #u_p[distance .< 5] .= NaN
    #v_p[distance .< 5] .= NaN
    #w_p[distance .< 5] .= NaN

    #ps = [Point3f(xi/8, yi/8, zi/2) for xi in x for yi in y for zi in z]
    #ns = [Vec3f(u_p[i], v_p[i], w_p[i]) for i in 1:length(u_p)]
    fig = Figure()
    ax = Axis(fig[1,1])
    #divergence[flags_tot .> 0] .= NaN
    hm = CairoMakie.heatmap!(ax, (x.-1), (z.-1).*4, divergence[:,31,:], colormap=:balance)
    arrows!(ax, (x.-1), (z.-1).*4, u_p[:,31,:], w_p[:,31,:].*4)
    Colorbar(fig[:, end+1], hm)
    #fig = arrows(ps, ns, fxaa=true)
    #save("$folder/quiver_rpca.png", fig)
    fig
end
main()
