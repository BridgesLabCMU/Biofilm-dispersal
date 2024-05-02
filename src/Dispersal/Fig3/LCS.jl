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
    folder = "/mnt/h/Dispersal/WT_replicate1_processed/Displacements/" 
    images_folder = dirname(dirname(folder))
    images = sort([f for f in readdir(images_folder, join=true) if occursin("isotropic.tif", f) && occursin("no_plank", f)], 
                  lt=natural)
    @show images
    img_size = size(TiffImages.load(images[1]; lazyio=true))
    files = sort([f for f in readdir(folder, join=true) if occursin("piv", f)], lt=natural)
    ntimepoints = length(files)
    x, y, z, u_dummy = load(files[1], "x", "y", "z", "u")
    height, width, nslices = size(u_dummy)
    u_tot = zeros(width, height, nslices, ntimepoints)
    v_tot = zeros(width, height, nslices, ntimepoints)
    w_tot = zeros(width, height, nslices, ntimepoints)
    for i in 1:length(files)
        u, v, w, flags = load(files[i], "u", "v", "w", "flags")
        u[flags .!= 0] .= 0
        v[flags .!= 0] .= 0
        w[flags .!= 0] .= 0
        @views u_tot[:,:,:,i] = permutedims(u[end:-1:1,:,:], (2,1,3))
        @views v_tot[:,:,:,i] = permutedims(v[end:-1:1,:,:], (2,1,3))
        @views w_tot[:,:,:,i] = permutedims(w[end:-1:1,:,:], (2,1,3))
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
    dx = 4 
    x0 = 1:dx:img_size[2]
    y0 = 1:dx:img_size[1]
    z0 = 1:dx:img_size[3]
    xlen = length(x0)
    ylen = length(y0)
    zlen = length(z0)
    dt = 0.1  
    dt2 = 1
    Tin = 5 
    Nsim = ntimepoints
    T = Nsim*dt2 
    yIC = zeros(Float32, 3, xlen, ylen, zlen) 
    solfor = zeros(Float64, xlen, ylen, zlen) 
    solbac = zeros(Float64, xlen, ylen, zlen)
    
    for i in 1:xlen
        for j in 1:ylen
            for k in 1:zlen
                @views yIC[:,i,j,k] = [x0[i], y0[j], z0[k]]
            end
        end
    end

    for m in 0:dt2:T
		@show m
        yin_for = yIC 
        for τ in m:dt:Tin+m
            yout = rk4singlestep((t,y) -> doublegyreVEC(t,y), dt, τ, yin_for)
            yin_for = yout
        end
        
		@views xT = yin_for[1,:,:,:]
        @views yT = yin_for[2,:,:,:]
        @views zT = yin_for[3,:,:,:]

        D = zeros(Float64, 3, 3)
        for i in 1:xlen
            for j in 1:ylen
                for k in 1:zlen
					@views dxTdx0 = i == 1 ? (xT[i+1,j,k] - xT[i,j,k])/dx :
							 i == xlen ? (xT[i,j,k] - xT[i-1,j,k])/dx :
							 (xT[i+1,j,k] - xT[i-1,j,k])/(2*dx)

					@views dxTdy0 = j == 1 ? (xT[i,j+1,k] - xT[i,j,k])/dx :
							 j == ylen ? (xT[i,j,k] - xT[i,j-1,k])/dx :
							 (xT[i,j+1,k] - xT[i,j-1,k])/(2*dx)

					@views dxTdz0 = k == 1 ? (xT[i,j,k+1] - xT[i,j,k])/dx :
							 k == zlen ? (xT[i,j,k] - xT[i,j,k-1])/dx :
							 (xT[i,j,k+1] - xT[i,j,k-1])/(2*dx)

					@views dyTdx0 = i == 1 ? (yT[i+1,j,k] - yT[i,j,k])/dx :
							 i == xlen ? (yT[i,j,k] - yT[i-1,j,k])/dx :
							 (yT[i+1,j,k] - yT[i-1,j,k])/(2*dx)

					@views dyTdy0 = j == 1 ? (yT[i,j+1,k] - yT[i,j,k])/dx :
							 j == ylen ? (yT[i,j,k] - yT[i,j-1,k])/dx :
							 (yT[i,j+1,k] - yT[i,j-1,k])/(2*dx)

					@views dyTdz0 = k == 1 ? (yT[i,j,k+1] - yT[i,j,k])/dx :
							 k == zlen ? (yT[i,j,k] - yT[i,j,k-1])/dx :
							 (yT[i,j,k+1] - yT[i,j,k-1])/(2*dx)

					@views dzTdx0 = i == 1 ? (zT[i+1,j,k] - zT[i,j,k])/dx :
							 i == xlen ? (zT[i,j,k] - zT[i-1,j,k])/dx :
							 (zT[i+1,j,k] - zT[i-1,j,k])/(2*dx)

					@views dzTdy0 = j == 1 ? (zT[i,j+1,k] - zT[i,j,k])/dx :
							 j == ylen ? (zT[i,j,k] - zT[i,j-1,k])/dx :
							 (zT[i,j+1,k] - zT[i,j-1,k])/(2*dx)

					@views dzTdz0 = k == 1 ? (zT[i,j,k+1] - zT[i,j,k])/dx :
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
                    solfor[i,j,k] = abs(1/Tin) * log(sqrt(maximum(eigvals(D'*D))))
                end
            end
        end
        solformin = minimum(solfor)
        solformax = maximum(solfor)
        solfor .-= solformin 
        solfor ./= solformax - solformin 

        yin_bac = yIC 
        for τ in m:-dt:-Tin+m
            yout = rk4singlestep((t,y) -> doublegyreVEC(t,y), -dt, τ, yin_bac)
            yin_bac = yout
        end
        @views xT = yin_bac[1,:,:,:]
        @views yT = yin_bac[2,:,:,:]
        @views zT = yin_bac[3,:,:,:]

        D = zeros(Float64, 3, 3)
        for i in 1:xlen
            for j in 1:ylen
                for k in 1:zlen
					@views dxTdx0 = i == 1 ? (xT[i+1,j,k] - xT[i,j,k])/dx :
							 i == xlen ? (xT[i,j,k] - xT[i-1,j,k])/dx :
							 (xT[i+1,j,k] - xT[i-1,j,k])/(2*dx)

					@views dxTdy0 = j == 1 ? (xT[i,j+1,k] - xT[i,j,k])/dx :
							 j == ylen ? (xT[i,j,k] - xT[i,j-1,k])/dx :
							 (xT[i,j+1,k] - xT[i,j-1,k])/(2*dx)

					@views dxTdz0 = k == 1 ? (xT[i,j,k+1] - xT[i,j,k])/dx :
							 k == zlen ? (xT[i,j,k] - xT[i,j,k-1])/dx :
							 (xT[i,j,k+1] - xT[i,j,k-1])/(2*dx)

					@views dyTdx0 = i == 1 ? (yT[i+1,j,k] - yT[i,j,k])/dx :
							 i == xlen ? (yT[i,j,k] - yT[i-1,j,k])/dx :
							 (yT[i+1,j,k] - yT[i-1,j,k])/(2*dx)

					@views dyTdy0 = j == 1 ? (yT[i,j+1,k] - yT[i,j,k])/dx :
							 j == ylen ? (yT[i,j,k] - yT[i,j-1,k])/dx :
							 (yT[i,j+1,k] - yT[i,j-1,k])/(2*dx)

					@views dyTdz0 = k == 1 ? (yT[i,j,k+1] - yT[i,j,k])/dx :
							 k == zlen ? (yT[i,j,k] - yT[i,j,k-1])/dx :
							 (yT[i,j,k+1] - yT[i,j,k-1])/(2*dx)

					@views dzTdx0 = i == 1 ? (zT[i+1,j,k] - zT[i,j,k])/dx :
							 i == xlen ? (zT[i,j,k] - zT[i-1,j,k])/dx :
							 (zT[i+1,j,k] - zT[i-1,j,k])/(2*dx)

					@views dzTdy0 = j == 1 ? (zT[i,j+1,k] - zT[i,j,k])/dx :
							 j == ylen ? (zT[i,j,k] - zT[i,j-1,k])/dx :
							 (zT[i,j+1,k] - zT[i,j-1,k])/(2*dx)

					@views dzTdz0 = k == 1 ? (zT[i,j,k+1] - zT[i,j,k])/dx :
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
                    solbac[i,j,k] = abs(1/Tin) * log(sqrt(maximum(eigvals(D'*D))))
                end
            end
        end
        solbacmin = minimum(solbac)
        solbacmax = maximum(solbac)
        solbac .-= solbacmin 
        solbac ./= solbacmax - solbacmin 
        save(folder*"FTLE_$(m).jld2", Dict("forward_LCS" => solfor, "backward_LCS" => solbac))
        solbac .= 0
        solfor .= 0
    end
end
main()
