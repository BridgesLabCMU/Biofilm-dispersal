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
                              (x, y, z, StepRangeLen(0, 1, ntimepoints))), 0)
    v_int = extrapolate(scale(interpolate(v_tot, BSpline(Cubic(Line(OnGrid())))), 
                              (x, y, z, StepRangeLen(0, 1, ntimepoints))), 0)
    w_int = extrapolate(scale(interpolate(w_tot, BSpline(Cubic(Line(OnGrid())))), 
                              (x, y, z, StepRangeLen(0, 1, ntimepoints))), 0)

    function doublegyreVEC(t, yin)
        @views shape = size(yin[1,:,:,:])
        @views x = vec(yin[1,:,:,:])
        @views y = vec(yin[2,:,:,:])
        @views z = vec(yin[3,:,:,:])
        t = zeros(Float64, size(x)) .+ t
        u = reshape(u_int.(x, y, z, t), shape)
        v = reshape(v_int.(x, y, z, t), shape)
        w = reshape(w_int.(x, y, z, t), shape)
        return permutedims(cat(u[:,:,:,:], v[:,:,:,:], w[:,:,:,:]; dims=4), (4,1,2,3))
    end

    # Constants
    Nsim = ntimepoints  
    dt = 0.1 
    T = Nsim

    # Set up a grid of particles
    dx = 5 # particle every voxel  
    x0 = 1:dx:img_size[2]
    y0 = 1:dx:img_size[1]
    z0 = 1:dx:img_size[3]
    xlen = length(x0)
    ylen = length(y0)
    zlen = length(z0)
    sol = zeros(Float32, 3, xlen, ylen, zlen) 

    for i in 1:xlen
        for j in 1:ylen
            for k in 1:zlen
                sol[:,i,j,k] = [x0[i], y0[j], z0[k]]
            end
        end
    end
    save(folder*"trajectories_0.jld2", Dict("trajectories" => sol))

    @inbounds for m in 0:T/dt
        time = m*dt
        @show time
        yout = rk4singlestep((t,y) -> doublegyreVEC(t,y), dt, time, sol)
        sol = yout
        save(folder*"particle_state_$(Int(m+1)).jld2", Dict("trajectories" => sol))
    end
end
main()
