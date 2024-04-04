# Load in u, v, w so they form 4D arrays (x, y, z, time)
#
# Interpolate the velocity fields in time and space
# Define a function (doublegyreVEC) which takes t, yin as arguments
# and returns dydt (i.e., the vector field) at yin and time t (based on interpolation)
# Define rk4singlestep

using FileIO
using Interpolations
using LinearAlgebra
using NaturalSort

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
    Nsim = 50
    dx = 0.7
    xvec = 
    yvec = 
    yIC = 
    dt=  
    dt2 = 
    Tin =
    T = Nsim*dt2

    solfor = []
    solbac = []

    for m in 0:dt2:T
        yin_for = yIC
        for i in m:dt:Tin
            yout_for = rk4singlestep(doublegyreVEC(t, y), dt, i, yin_for)
            yin_for = yout_for
        end

        xT = yin_for[1]
        yT = yin_for[2]

        dxTdx0, dxTdy0 = gradient(xT, dx, dx)
        dyTdx0, dyTdy0 = gradient(yT, dx, dx)

        D = 
        sigma = zeros(Float64, size(xT))
        for i in 1:length(xvec)
            for j in 1:length(yvec)
                D[0,0] = dxTdx0[j,i]
                D[0,1] = dxTdy0[j,i]
                D[1,0] = dyTdx0[j,i]
                D[1,1] = dyTdy0[j,i]
                sigma[j,i] = abs(1/Tin) * max(eigvals(D'*D))
            end
        end
        sigma = (sigma .- minimum(sigma)) ./ (maximum(sigma) - minimum(sigma))
        push!(solfor, sigma)

        yin_bac = yIC

        for 
    end
end
