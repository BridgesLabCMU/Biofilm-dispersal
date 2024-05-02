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

function calculate_radial_component(u, v, w, dx,dy,dz, central_point)
    nx, ny, nz = size(u)
    x = reshape(1:nx, nx, 1, 1)
    y = reshape(1:ny, 1, ny, 1)
    z = reshape(1:nz, 1, 1, nz)
    dx = x .- central_point[1]
    dy = y .- central_point[2]
    dz = z .- central_point[3]
    magnitudes = sqrt.(dx.^2 + dy.^2 + dz.^2)
    dx_norm = dx ./ magnitudes
    dy_norm = dy ./ magnitudes
    dz_norm = dz ./ magnitudes
    radial_components = u .* dx_norm + v .* dy_norm + w .* dz_norm
    return radial_components
end

function compute_com(mask, vector_field)
	mask = label_components(mask)
	areas = component_lengths(mask)
	max_label = argmax(areas[1:end])
	mask[mask .!= max_label] .= 0
	mask[mask .== max_label] .= 1 
    center_of_mass = component_centroids(mask)[1]
    relative_com = [x for x in center_of_mass] ./ [x for x in size(mask)] 
    center_of_mass = permutedims(relative_com, (2,1,3)) .* [x for x in size(vector_field)] .- 1
    center_of_mass[3] = 0
	return center_of_mass
end

function main()
    folder = "/mnt/h/Dispersal/WT_replicate1_processed/Displacements/"
    images_folder = "/mnt/h/Dispersal/WT_replicate1_processed/"
    mask_files = sort([f for f in readdir(images_folder, join=true) if occursin("mask_isotropic", f)], lt=natural)
    mask_t1 = load(mask_files[1])

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
    center_of_mass = compute_com(mask_t1, u_p)

    u_p[flags_tot .> 0] .= NaN 
    v_p[flags_tot .> 0] .= NaN  
    w_p[flags_tot .> 0] .= NaN
    radial_component = calculate_radial_component(u_p, v_p, w_p.*4, 8, 8, 8, center_of_mass)

    fig = Figure()
    ax = Axis(fig[1,1])
    hm = CairoMakie.heatmap!(ax, (x.-1), (z.-1).*4, radial_component[:,31,:], colormap=:balance)
    arrows!(ax, (x.-1), (z.-1).*4, u_p[:,31,:], w_p[:,31,:].*4)
    Colorbar(fig[:, end+1], hm)
    fig
end
main()
