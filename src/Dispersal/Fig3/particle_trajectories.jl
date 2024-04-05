using FileIO
using Interpolations
using LinearAlgebra
using NaturalSort
using TiffImages

function rk4singlestep(f, dt, t, y)
    k1 = f(t, y)
    k2 = f(t + dt/2, y + dt/2*k1)
    k3 = f(t + dt/2, y + dt/2*k2)
    k4 = f(t + dt, y + dt*k3)
    return y + dt/6*(k1 + 2*k2 + 2*k3 + k4)
end

function main()
    # Load in the data
    folder = "/mnt/h/Dispersal/WT_replicate1_processed/Displacements/" 
    images_folder = dirname(dirname(folder))
    images = sort([f for f in readdir(images_folder, join=true) if occursin("isotropic.tif", f) && occursin("noplank", f)], 
                  lt=natural)
    img_size = size(load(images[1]); lazyio=true)
    files = sort([f for f in readdir(folder, join=true) if occursin("piv", f)], lt=natural)
    ntimepoints = length(files)
    x, y, z, u_dummy = load(files[1], "x", "y", "z", "u")
    height, width, nslices = size(u_dummy)
    u_tot = zeros(width, height, nslices, ntimepoints)
    v_tot = zeros(width, height, nslices, ntimepoints)
    w_tot = zeros(width, height, nslices, ntimepoints)
    for i in 1:length(files)
        u, v, w = load(files[i], "u", "v", "w")
        u_tot[:,:,:,i] = permutedims(u, (2,1,3))
        v_tot[:,:,:,i] = permutedims(v, (2,1,3))
        w_tot[:,:,:,i] = permutedims(w, (2,1,3))
    end

    w_tot .*= 4
    z .-= 1
    z .*= 4
    z .+= 1

    u_int = extrapolate(interpolate((x, y, z, 1:ntimepoints), u_tot, Gridded(Linear())), 0)
    v_int = extrapolate(interpolate((x, y, z, 1:ntimepoints), v_tot, Gridded(Linear())), 0)
    w_int = extrapolate(interpolate((x, y, z, 1:ntimepoints), w_tot, Gridded(Linear())), 0)

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
    dx = 5 
    x0 = 1:dx:img_size[2]
    y0 = 1:dx:img_size[1]
    z0 = 1:dx:img_size[3]
    xlen = length(x0)
    ylen = length(y0)
    zlen = length(z0)
    xT = copy(x0) 
    yT = copy(y0) 
    zT = copy(z0) 
    dt = 0.1  
    dt2 = Delta
    Tin = 40
    T = Nsim*dt2
    sol_t = zeros(Float64, xlen, ylen, zlen, Nsim) 

    for i in 1:xlen
        for j in 1:ylen
            for k in 1:zlen
                sol_t[i,j,k,1] = [x0[i], y0[j], z0[k]]
            end
        end
    end

    for m in 2:T
        for i in 1:xlen
            for j in 1:ylen 
                for k in 1:zlen
                    for τ in m-1:dt:Tin
                        sol = rk4singlestep((t,y) -> doublegyreVEC(t, y), dt, 
                                            τ, [x0[i],y0[j],z0[k]])
                    end
                    xT[i,j,k] = sol[1]
                    yT[i,j,k] = sol[2]
                    zT[i,j,k] = sol[3]
                    sol_t[i,j,k,m] = [sol[1],sol[2],sol[3]]
                end
            end
        end
    end
    save(folder*"trajectories.jld2", Dict("trajectories" => sol_t))
end
