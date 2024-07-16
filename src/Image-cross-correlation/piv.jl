module PIV

export piv

using FileIO
using TiffImages
using Interpolations
using NaturalSort
using FFTW
using StatsBase
using ImageCore: channelview
using Images: imresize
using ImageFiltering: mapwindow, Fill
using NaNStatistics
using Memoize
using Revise

include("validation.jl")
include("crossCorr.jl")
include("utils.jl")

function multipass(x_old, y_old, z_old, u_old, v_old, w_old, s2n, volume1, volume2, i, params)
    # Compute new coordinates
    x_new, y_new, z_new = compute_coords(size(volume1), params.windows[i], params.overlaps[i])
    
    if x_old != nothing
        itp_u = interpolate(u_old, BSpline(Cubic(Line(OnGrid()))))
        itp_v = interpolate(v_old, BSpline(Cubic(Line(OnGrid()))))
        itp_w = interpolate(w_old, BSpline(Cubic(Line(OnGrid()))))
        sitp_u = extrapolate(scale(itp_u, (y_old, x_old, z_old)), Line())
        sitp_v = extrapolate(scale(itp_v, (y_old, x_old, z_old)), Line())
        sitp_w = extrapolate(scale(itp_w, (y_old, x_old, z_old)), Line())
        u_old = sitp_u(y_new, x_new, z_new) 
        v_old = sitp_v(y_new, x_new, z_new) 
        w_old = sitp_w(y_new, x_new, z_new) 

        # Create deformation field 
        h, w, d = size(volume1)
        xt = 1:w 
        yt = 1:h
        zt = 1:d
        itp_u_img = extrapolate(scale(interpolate(u_old, BSpline(Cubic(Line(OnGrid())))), 
                                      (y_new, x_new, z_new)), Line())
        itp_v_img = extrapolate(scale(interpolate(v_old, BSpline(Cubic(Line(OnGrid())))), 
                                      (y_new, x_new, z_new)), Line())
        itp_w_img = extrapolate(scale(interpolate(w_old, BSpline(Cubic(Line(OnGrid())))), 
                                      (y_new, x_new, z_new)), Line())
        ut = itp_u_img(yt, xt, zt) 
        vt = itp_v_img(yt, xt, zt) 
        wt = itp_w_img(yt, xt, zt)

        # Perform symmetric deformation
        itp_vol1 = extrapolate(interpolate(volume1, BSpline(Cubic(Line(OnGrid())))), Line())
        itp_vol2 = extrapolate(interpolate(volume2, BSpline(Cubic(Line(OnGrid())))), Line())
        volume1 = itp_vol1.(yt .- vt ./ 2, reshape(xt,1,:,1) .- ut ./ 2, 
                           reshape(zt,1,1,:) .- wt ./ 2)
        volume2 = itp_vol2.(yt .+ vt ./ 2, reshape(xt,1,:,1) .+ ut ./ 2, 
                           reshape(zt,1, 1,:) .+ wt ./ 2)
    end

    shape = field_shape(size(volume1), params.windows[i], params.overlaps[i])
    # Pre-allocate
    u_new = zeros(Float32, shape)
    v_new = zeros(Float32, shape)
    w_new = zeros(Float32, shape)
    s2n = zeros(Float32, shape)
    compute_displacements!(volume1, volume2, params.windows[i], x_new, y_new, z_new, 
                           params.overlaps[i], params.s2n_method,
                           u_new, v_new, w_new, s2n)
    if u_old != nothing
        u_new += u_old
        v_new += v_old
        w_new += w_old
    end
    return x_new, y_new, z_new, u_new, v_new, w_new, s2n
end

function piv(params)
    
    # Get an array of timelapse files
    files = sort([f for f in readdir(params.image_folder, join=true) if occursin("noplank", f) && !occursin("isotropic", f)],
                 lt=natural)

    # Create a folder to save the results if it doesn't exist
    if !isdir("Displacements")
        mkdir("Displacements")
    end

    for t in 2:length(files)
        volume1 = imresize(Float32.(channelview(TiffImages.load(files[t-1]))), ratio=(1,1,(0.3/0.065)/4))
        volume2 = imresize(Float32.(channelview(TiffImages.load(files[t]))), ratio=(1,1,(0.3/0.065)/4))
        x = nothing
        y = nothing
        z = nothing
        u = nothing
        v = nothing
        w = nothing
        s2n = nothing
        flags = nothing

        for i in 1:params.mpass
            x, y, z, u, v, w, s2n = multipass(x, y, z, u, v, w, s2n, 
                                              volume1, volume2, i, params)
            flags = zeros(Bool, size(u))
            validate!(flags, u, v, w, s2n, i, params)
            u, v, w = replace_outliers(u, v, w, flags, params; 
                                       max_iter=params.max_iter, 
                                       tol=params.tol,
                                       kernel_size=params.kernel_size)
        end
        # Transform the coordinates
        x, y, z, u, v, w = transform_coords(x, y, z, u, v, w)
        # Save
        save(params.image_folder*"/Displacements/piv_results_$(t-1).jld2", 
             Dict("x" => x, "y" => y, "z" => z, "u" => u, 
                  "v" => v, "w" => w, "s2n" => s2n, "flags" => flags))
    end
end

end # module
