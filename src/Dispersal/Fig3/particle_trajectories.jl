using FileIO
using Interpolations
using LinearAlgebra
using NaturalSort
using TiffImages

function rk4singlestep(f, dt, t, y)
    k1 = f(t, y)
    k2 = f(t + dt/2, y + k1.*dt/2)
    k3 = f(t + dt/2, y + k2.*dt/2)
    k4 = f(t + dt, y + k3.*dt)
    return y .+ dt/6 .* (k1 + k2.*2 + k3.*2 + k4)
end

function main()
    # Load in the data
    folder = "/mnt/h/Dispersal/WT_replicate1_processed/Displacements/" 
    images_folder = dirname(dirname(folder))
    images = sort([f for f in readdir(images_folder, join=true) if occursin("isotropic.tif", f) && occursin("no_plank", f)], 
                  lt=natural)
    img_size = size(TiffImages.load(images[1]; lazyio=true))
    files = sort([f for f in readdir(folder, join=true) if occursin("piv", f)], lt=natural)
    ntimepoints = length(files)
    x, y, z, u_dummy = load(files[1], "x", "y", "z", "u")
    height, width, nslices = size(u_dummy)
    u_tot = zeros(width, height, nslices, ntimepoints)
    v_tot = zeros(width, height, nslices, ntimepoints)
    w_tot = zeros(width, height, nslices, ntimepoints)
    for i in 1:length(files)
        u, v, w = load(files[i], "u", "v", "w")
        u_tot[:,:,:,i] = permutedims(u[end:-1:1,:,:], (2,1,3))
        v_tot[:,:,:,i] = permutedims(v[end:-1:1,:,:], (2,1,3))
        w_tot[:,:,:,i] = permutedims(w[end:-1:1,:,:], (2,1,3))
    end

    w_tot .*= 4
    z = z .- 1
    z = z .* 4
    z = z .+ 1
    y = y[end:-1:1]
    
    u_int = extrapolate(scale(interpolate(u_tot, BSpline(Cubic(Line(OnGrid())))), 
                              (x, y, z, StepRangeLen(1, 1, ntimepoints))), Line())
    v_int = extrapolate(scale(interpolate(v_tot, BSpline(Cubic(Line(OnGrid())))), 
                              (x, y, z, StepRangeLen(1, 1, ntimepoints))), Line())
    w_int = extrapolate(scale(interpolate(w_tot, BSpline(Cubic(Line(OnGrid())))), 
                              (x, y, z, StepRangeLen(1, 1, ntimepoints))), Line())

    function doublegyreVEC(t, yin)
        x, y, z = yin
        u = u_int[x, y, z, t]
        v = v_int[x, y, z, t]
        w = w_int[x, y, z, t]
        return [u, v, w]
    end

    # Constants
    Delta = 1
    Nsim = ntimepoints  
    dx = 30 
    x0 = 1:dx:img_size[2]
    y0 = 1:dx:img_size[1]
    z0 = 1:dx:img_size[3]
    xlen = length(x0)
    ylen = length(y0)
    zlen = length(z0)
    dt = 0.1  
    dt2 = Delta
    Tin = 40
    T = Nsim*dt2
    sol_t = zeros(Float64, 3, xlen, ylen, zlen, Nsim) 

    for i in 1:xlen
        for j in 1:ylen
            for k in 1:zlen
                sol_t[1:3, i,j,k,1] = [x0[i], y0[j], z0[k]]
            end
        end
    end

    @inbounds for m in 2:T
        @show m
        @inbounds for i in 1:xlen
            @inbounds for j in 1:ylen 
                @inbounds for k in 1:zlen
                    @inbounds for τ in m-1:dt:Tin+m
                        sol_t[:,i,j,k,m] = rk4singlestep((t,y) -> doublegyreVEC(t, y), dt, 
                                            τ, sol_t[:,i,j,k,m])
                    end
                end
            end
        end
    end
    save(folder*"trajectories.jld2", Dict("trajectories" => sol_t))
end
main()
