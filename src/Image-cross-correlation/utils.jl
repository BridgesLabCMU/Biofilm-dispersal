mutable struct PIVParams
    image_folder::String
    mpass::Int
    windows::Array{Tuple{Int, Int, Int}}
    overlaps::Array{Tuple{Int, Int, Int}}
    u_thresh::Tuple{Float64, Float64}
    v_thresh::Tuple{Float64, Float64}
    w_thresh::Tuple{Float64, Float64}
    sig2noise_thresh::Float64
    median_thresh::Int
    median_size::Int
    std_thresh::Int
    std_size::Int
    s2n_thresh::Float64
    s2n_method::String
    replace_method::String
    max_iter::Int
    tol::Float64
    kernel_size::Int
    median_validate::Bool
    std_validate::Bool
    global_validate::Bool
    s2n_validate::Bool
end

function pivSettings(;image_folder="", mpass=3, windows=[(64,64,16), (32,32,8), (16,16,4)], 
                        overlaps=[(32,32,8), (16,16,4), (8,8,2)], 
                        u_thresh=(-32.0, 32.0), v_thresh=(-32.0, 32.0), 
                        w_thresh=(-8.0, 8.0), 
                        sig2noise_thresh=1.0, median_thresh=3, 
                        median_size=1, std_thresh=10, std_size=1, s2n_thresh=1.0,
                        s2n_method="peak2mean", replace_method="localmean",  
                        max_iter=9, tol=1e-3, kernel_size=2, median_validate=true, 
                        std_validate=true, global_validate=true, s2n_validate=true)
    return PIVParams(image_folder, mpass, windows, overlaps, u_thresh, v_thresh, w_thresh, 
                     sig2noise_thresh, median_thresh, median_size, std_thresh,
                     std_size, s2n_thresh, s2n_method, replace_method, max_iter, tol,
                     kernel_size, median_validate, std_validate, global_validate, s2n_validate)
end

function transform_coords(x, y, z, u, v, w)
    y = y[end:-1:1]
    v .*= -1
    return x, y, z, u, v, w
end

function field_shape(image_shape, window_size, overlap)
    return tuple([(i - w) ÷ (w - o) + 1 for 
                  (i, w, o) in zip(image_shape, window_size, overlap)]...)
end

function compute_coords(image_shape, window_size, overlap)
    shape = field_shape(image_shape, window_size, overlap)
    x = (0:shape[2]-1) .* (window_size[2] - overlap[2]) .+ window_size[2] ÷ 2 .+ 1
    y = (0:shape[1]-1) .* (window_size[1] - overlap[1]) .+ window_size[1] ÷ 2 .+ 1
    z = (0:shape[3]-1) .* (window_size[3] - overlap[3]) .+ window_size[3] ÷ 2 .+ 1
    return x, y, z
end
