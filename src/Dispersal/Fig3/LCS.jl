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
                              (x, y, z, StepRangeLen(1, 1, ntimepoints))), 0)
    v_int = extrapolate(scale(interpolate(v_tot, BSpline(Cubic(Line(OnGrid())))), 
                              (x, y, z, StepRangeLen(1, 1, ntimepoints))), 0)
    w_int = extrapolate(scale(interpolate(w_tot, BSpline(Cubic(Line(OnGrid())))), 
                              (x, y, z, StepRangeLen(1, 1, ntimepoints))), 0)

    function doublegyreVEC(t, yin)
        x, y, z = yin
        u = u_int[x, y, z, t]
        v = v_int[x, y, z, t]
        w = w_int[x, y, z, t]
        return [u, v, w]
    end

    # Constants
    dx = 1
    x0 = 1:dx:img_size[2]
    y0 = 1:dx:img_size[1]
    z0 = 1:dx:img_size[3]
    xlen = length(x0)
    ylen = length(y0)
    zlen = length(z0)
    xT = zeros(Float64, xlen, ylen, zlen)
    yT = zeros(Float64, xlen, ylen, zlen)
    zT = zeros(Float64, xlen, ylen, zlen)
    dt = 0.1  
    Tin = 40
    T = ntimepoints 

    solfor = zeros(Float64, xlen, ylen, zlen, T) 
    solbac = zeros(Float64, xlen, ylen, zlen, T)

    for m in 1:T
        @show m
        for i in 1:xlen
            for j in 1:ylen 
                for k in 1:zlen
                    for τ in m:dt:Tin+m
                        sol = rk4singlestep((t,y) -> doublegyreVEC(t, y), dt, 
                                            τ, [x0[i],y0[j],z0[k]])
                    end
                    xT[i,j,k] = sol[1]
                    yT[i,j,k] = sol[2]
                    zT[i,j,k] = sol[3]
                end
            end
        end

        D = zeros(Float64, 3, 3)
        for i in 1:xlen
            for j in 1:ylen
                for k in 1:zlen
					dxTdx0 = i == 1 ? (xT[i+1,j,k] - xT[i,j,k])/dx :
							 i == xlen ? (xT[i,j,k] - xT[i-1,j,k])/dx :
							 (xT[i+1,j,k] - xT[i-1,j,k])/(2*dx)

					dxTdy0 = j == 1 ? (xT[i,j+1,k] - xT[i,j,k])/dx :
							 j == ylen ? (xT[i,j,k] - xT[i,j-1,k])/dx :
							 (xT[i,j+1,k] - xT[i,j-1,k])/(2*dx)

					dxTdz0 = k == 1 ? (xT[i,j,k+1] - xT[i,j,k])/dx :
							 k == zlen ? (xT[i,j,k] - xT[i,j,k-1])/dx :
							 (xT[i,j,k+1] - xT[i,j,k-1])/(2*dx)

					dyTdx0 = i == 1 ? (yT[i+1,j,k] - yT[i,j,k])/dx :
							 i == xlen ? (yT[i,j,k] - yT[i-1,j,k])/dx :
							 (yT[i+1,j,k] - yT[i-1,j,k])/(2*dx)

					dyTdy0 = j == 1 ? (yT[i,j+1,k] - yT[i,j,k])/dx :
							 j == ylen ? (yT[i,j,k] - yT[i,j-1,k])/dx :
							 (yT[i,j+1,k] - yT[i,j-1,k])/(2*dx)

					dyTdz0 = k == 1 ? (yT[i,j,k+1] - yT[i,j,k])/dx :
							 k == zlen ? (yT[i,j,k] - yT[i,j,k-1])/dx :
							 (yT[i,j,k+1] - yT[i,j,k-1])/(2*dx)

					dzTdx0 = i == 1 ? (zT[i+1,j,k] - zT[i,j,k])/dx :
							 i == xlen ? (zT[i,j,k] - zT[i-1,j,k])/dx :
							 (zT[i+1,j,k] - zT[i-1,j,k])/(2*dx)

					dzTdy0 = j == 1 ? (zT[i,j+1,k] - zT[i,j,k])/dx :
							 j == ylen ? (zT[i,j,k] - zT[i,j-1,k])/dx :
							 (zT[i,j+1,k] - zT[i,j-1,k])/(2*dx)

					dzTdz0 = k == 1 ? (zT[i,j,k+1] - zT[i,j,k])/dx :
							 k == zlen ? (zT[i,j,k] - zT[i,j,k-1])/dx :
							 (zT[i,j,k+1] - zT[i,j,k-1])/(2*dx)
                    D[1,1] = dxTdx0
                    D[1,2] = dxTdy0
                    D[1,3] = dxTdz0
                    D[2,1] = dyTdx0
                    D[2,2] = dyTdy0
                    D[2,3] = dyTdz0
                    D[3,1] = dzTdx0
                    D[3,2] = dzTdy0
                    D[3,3] = dzTdz0
                    solbac[i,j,k,m] = abs(1/Tin) * maximum(eigvals(D'*D))
                end
            end
        end
        @views solformin = minimum(solfor[:,:,:,m])
        @views solformax = maximum(solfor[:,:,:,m])
        solfor[:,:,:,m] .-= solformin 
        solfor[:,:,:,m] ./= solformax - solformin 

        for i in 1:xlen
            for j in 1:ylen
                for k in 1:zlen
                    for τ in m:-dt:-Tin+m
                        sol = rk4singlestep((t,y) -> doublegyreVEC(t, y), -dt, 
                                            τ, [x0[i],y0[j],z0[k]])
                    end
                    xT[i,j,k] = sol[1]
                    yT[i,j,k] = sol[2]
                    zT[i,j,k] = sol[3]
                end
            end
        end

        D = zeros(Float64, 3, 3)
        for i in 1:xlen
            for j in 1:ylen
                for k in 1:zlen
					dxTdx0 = i == 1 ? (xT[i+1,j,k] - xT[i,j,k])/dx :
							 i == xlen ? (xT[i,j,k] - xT[i-1,j,k])/dx :
							 (xT[i+1,j,k] - xT[i-1,j,k])/(2*dx)

					dxTdy0 = j == 1 ? (xT[i,j+1,k] - xT[i,j,k])/dx :
							 j == ylen ? (xT[i,j,k] - xT[i,j-1,k])/dx :
							 (xT[i,j+1,k] - xT[i,j-1,k])/(2*dx)

					dxTdz0 = k == 1 ? (xT[i,j,k+1] - xT[i,j,k])/dx :
							 k == zlen ? (xT[i,j,k] - xT[i,j,k-1])/dx :
							 (xT[i,j,k+1] - xT[i,j,k-1])/(2*dx)

					dyTdx0 = i == 1 ? (yT[i+1,j,k] - yT[i,j,k])/dx :
							 i == xlen ? (yT[i,j,k] - yT[i-1,j,k])/dx :
							 (yT[i+1,j,k] - yT[i-1,j,k])/(2*dx)

					dyTdy0 = j == 1 ? (yT[i,j+1,k] - yT[i,j,k])/dx :
							 j == ylen ? (yT[i,j,k] - yT[i,j-1,k])/dx :
							 (yT[i,j+1,k] - yT[i,j-1,k])/(2*dx)

					dyTdz0 = k == 1 ? (yT[i,j,k+1] - yT[i,j,k])/dx :
							 k == zlen ? (yT[i,j,k] - yT[i,j,k-1])/dx :
							 (yT[i,j,k+1] - yT[i,j,k-1])/(2*dx)

					dzTdx0 = i == 1 ? (zT[i+1,j,k] - zT[i,j,k])/dx :
							 i == xlen ? (zT[i,j,k] - zT[i-1,j,k])/dx :
							 (zT[i+1,j,k] - zT[i-1,j,k])/(2*dx)

					dzTdy0 = j == 1 ? (zT[i,j+1,k] - zT[i,j,k])/dx :
							 j == ylen ? (zT[i,j,k] - zT[i,j-1,k])/dx :
							 (zT[i,j+1,k] - zT[i,j-1,k])/(2*dx)

					dzTdz0 = k == 1 ? (zT[i,j,k+1] - zT[i,j,k])/dx :
							 k == zlen ? (zT[i,j,k] - zT[i,j,k-1])/dx :
							 (zT[i,j,k+1] - zT[i,j,k-1])/(2*dx)
                    D[1,1] = dxTdx0
                    D[1,2] = dxTdy0
                    D[1,3] = dxTdz0
                    D[2,1] = dyTdx0
                    D[2,2] = dyTdy0
                    D[2,3] = dyTdz0
                    D[3,1] = dzTdx0
                    D[3,2] = dzTdy0
                    D[3,3] = dzTdz0
                    solbac[i,j,k,m] = abs(1/Tin) * maximum(eigvals(D'*D))
                end
            end
        end
        @views solbacmin = minimum(solbac[:,:,:,m])
        @views solbacmax = maximum(solbac[:,:,:,m])
        solbac[:,:,:,m] .-= solbacmin 
        solbac[:,:,:,m] ./= solbacmax - solbacmin 
    end
    save(folder*"FTLE.jld2", Dict("forward_LCS" => solfor, "backward_LCS" => solbac))
end
main()
