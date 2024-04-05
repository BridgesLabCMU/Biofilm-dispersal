"""
Procedure:
1. Load in vector fields as 4D arrays (x,y,z,t) 
2. Interpolate them to get continuous vector fields in space and time
3. Define a function that for any configuration of particles, returns the velocity at that point
4. Define constants and initial conditions
5. Starting with the initial configuration (x, y, z meshgrid), step it forward in time using RK4
6. Retrieve the x, y, z coords after the timestep
7. Calculate the spatial gradient of the x, y, z coords
8. For each x, y, z position, calculate its deformation tensor
9. Calculate the maximum singular value and the FTLE at that point and add it to the forward LCS


Currently the code is completely vectorized aside from the calculation of the FTLE at each point.
Is it possible to get rid of all vectorization?





"""

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
    dx = 1
    x0 = 1:dx:img_size[2]
    y0 = 1:dx:img_size[1]
    z0 = 1:dx:img_size[3]
    xT = zeros(Float64, length(x0))
    yT = zeros(Float64, length(y0))
    zT = zeros(Float64, length(z0))
    dt = 0.1  
    dt2 = Delta
    Tin = 40
    T = Nsim*dt2

    solfor = zeros(Float64, length(x0), length(y0), length(z0), Nsim) 
    solbac = zeros(Float64, length(x0), length(y0), length(z0), Nsim)

    for m in 1:dt2:T
        for i in 1:length(x0)
            for j in 1:length(y0) 
                for k in 1:length(z0)
                    for τ in m:dt:Tin+m
                        sol = rk4singlestep(doublegyreVEC(t, y), dt, 
                                            τ, [x0[i],y0[j],z0[k]])
                    end
                    xT[i] = sol[1]
                    yT[j] = sol[2]
                    zT[k] = sol[3]
                end
            end
        end

        D = zeros(Float64, 3, 3)
        for i in 1:length(x0)
            for j in 1:length(y0)
                for k in 1:length(z0)
                    # Calculate gradients here
                    D[1,1] = dxTdx0[j,i,k]
                    D[1,2] = dxTdy0[j,i,k]
                    D[1,3] = dxTdz0[j,i,k]
                    D[2,1] = dyTdx0[j,i,k]
                    D[2,2] = dyTdy0[j,i,k]
                    D[2,3] = dyTdz0[j,i,k]
                    D[3,1] = dzTdx0[j,i,k]
                    D[3,2] = dzTdy0[j,i,k]
                    D[3,3] = dzTdz0[j,i,k]
                    solfor[i,j,k,m] = abs(1/Tin) * max(eigvals(D'*D))
                end
            end
        end
        solfor[:,:,:,m] = (solfor[:,:,:,m] .- minimum(solfor[:,:,:,m])) ./ (maximum(solfor[:,:,:,m]) - minimum(solfor[:,:,:,m]))

        for i in 1:length(x0)
            for j in 1:length(y0) 
                for k in 1:length(z0)
                    for τ in m:-dt:-Tin+m
                        sol = rk4singlestep(doublegyreVEC(t, y), -dt, 
                                            τ, [x0[i],y0[j],z0[k]])
                    end
                    xT[i] = sol[1]
                    yT[j] = sol[2]
                    zT[k] = sol[3]
                end
            end
        end

        D = zeros(Float64, 3, 3)
        for i in 1:length(x0)
            for j in 1:length(y0)
                for k in 1:length(z0)
                    # Calculate gradients here
                    D[1,1] = dxTdx0[j,i,k]
                    D[1,2] = dxTdy0[j,i,k]
                    D[1,3] = dxTdz0[j,i,k]
                    D[2,1] = dyTdx0[j,i,k]
                    D[2,2] = dyTdy0[j,i,k]
                    D[2,3] = dyTdz0[j,i,k]
                    D[3,1] = dzTdx0[j,i,k]
                    D[3,2] = dzTdy0[j,i,k]
                    D[3,3] = dzTdz0[j,i,k]
                    solbac[i,j,k,m] = abs(1/Tin) * max(eigvals(D'*D))
                end
            end
        end
        solbac[:,:,:,m] = (solbac[:,:,:,m] .- minimum(solbac[:,:,:,m])) ./ (maximum(solbac[:,:,:,m]) - minimum(solbac[:,:,:,m]))
    end
    save(folder*"FTLE.jld2", Dict("solfor" => solfor, "solbac" => solbac))
end
