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
    dx = 1
    xvec = 1:dx:img_size[2]
    yvec = 1:dx:img_size[1]
    zvec = 1:dx:img_size[3]
    zvec .-= 1
    zvec .*= 4
    zvec .+= 1
    yIC = zeros(Float64, (3, length(yvec), length(xvec), length(zvec)))
    dt = 0.1  
    dt2 = Delta
    Tin = 40
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
        zT = yin_for[3]

        dxTdx0, dxTdy0, dxTdz0 = gradient(xT, dx, dx)
        dyTdx0, dyTdy0, dyTdz0 = gradient(yT, dx, dx)
        dzTdx0, dzTdy0, dzTdz0 = gradient(zT, dx, dx)

        D = zeros(Float64, 3, 3)
        sigma = zeros(Float64, size(xT))
        for i in 1:length(xvec)
            for j in 1:length(yvec)
                for k in 1:length(zvec)
                    D[1,1] = dxTdx0[j,i,k]
                    D[1,2] = dxTdy0[j,i,k]
                    D[1,3] = dxTdz0[j,i,k]
                    D[2,1] = dyTdx0[j,i,k]
                    D[2,2] = dyTdy0[j,i,k]
                    D[2,3] = dyTdz0[j,i,k]
                    D[3,1] = dzTdx0[j,i,k]
                    D[3,2] = dzTdy0[j,i,k]
                    D[3,3] = dzTdz0[j,i,k]
                    sigma[j,i,k] = abs(1/Tin) * max(eigvals(D'*D))
                end
            end
        end
        sigma = (sigma .- minimum(sigma)) ./ (maximum(sigma) - minimum(sigma))
        push!(solfor, sigma)

        yin_bac = yIC

        for i in m:-dt:-Tin+m
            yout = rk4singlestep(doublegyreVEC(t, y), -dt, i, yin_bac)
            yin_bac = yout
        end

        xT = yin_bac[1]
        yT = yin_bac[2]
        zT = yin_bac[3]

        dxTdx0, dxTdy0, dxTdz0 = gradient(xT, dx, dx)
        dyTdx0, dyTdy0, dyTdz0 = gradient(yT, dx, dx)
        dzTdx0, dzTdy0, dzTdz0 = gradient(zT, dx, dx)

        D = zeros(Float64, 3, 3)
        sigma = zeros(Float64, size(xT))
        for i in 1:length(xvec)
            for j in 1:length(yvec)
                for k in 1:length(zvec)
                    D[1,1] = dxTdx0[j,i,k]
                    D[1,2] = dxTdy0[j,i,k]
                    D[1,3] = dxTdz0[j,i,k]
                    D[2,1] = dyTdx0[j,i,k]
                    D[2,2] = dyTdy0[j,i,k]
                    D[2,3] = dyTdz0[j,i,k]
                    D[3,1] = dzTdx0[j,i,k]
                    D[3,2] = dzTdy0[j,i,k]
                    D[3,3] = dzTdz0[j,i,k]
                    sigma[j,i,k] = abs(1/Tin) * max(eigvals(D'*D))
                end
            end
        end
        sigma = (sigma .- minimum(sigma)) ./ (maximum(sigma) - minimum(sigma))
        push!(solbac, sigma)
    end
    save(folder*"solfor.jld2", Dict("solfor" => solfor, "solbac" => solbac))
end
