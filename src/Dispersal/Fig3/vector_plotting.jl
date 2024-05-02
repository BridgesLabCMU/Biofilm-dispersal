using FileIO 
using NPZ
using CairoMakie
using NaturalSort
using StatsBase
using NaNStatistics

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

    u_p[flags_tot .> 0] .= NaN 
    v_p[flags_tot .> 0] .= NaN  
    w_p[flags_tot .> 0] .= NaN
    divergence = calculate_divergence(u_p, v_p, w_p.*4, 8, 8, 8)

    fig = Figure()
    ax = Axis(fig[1,1])
    hm = CairoMakie.heatmap!(ax, (x.-1), (z.-1).*4, divergence[:,31,:], colormap=:balance)
    arrows!(ax, (x.-1), (z.-1).*4, u_p[:,31,:], w_p[:,31,:].*4)
    Colorbar(fig[:, end+1], hm)
    fig
end
main()
