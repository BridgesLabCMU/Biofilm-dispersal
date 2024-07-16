function validate!(flags, u, v, w, s2n, i, params)
    if params.median_validate
        um = mapwindow(nanmedian, u, ntuple(x->params.median_size*2+1, 3), 
                       Fill(NaN))
        vm = mapwindow(nanmedian, v, ntuple(x->params.median_size*2+1, 3), 
                       Fill(NaN))
        wm = mapwindow(nanmedian, w, ntuple(x->params.median_size*2+1, 3),
                       Fill(NaN))
        flags[abs.(u .- um) .> params.median_thresh] .= true
        flags[abs.(v .- vm) .> params.median_thresh] .= true
        flags[abs.(w .- wm) .> params.median_thresh] .= true
    end
    if params.std_validate
        ustd = mapwindow(std, u, ntuple(x->params.std_size*2+1, 3))
        vstd = mapwindow(std, v, ntuple(x->params.std_size*2+1, 3))
        wstd = mapwindow(std, w, ntuple(x->params.std_size*2+1, 3))
        flags[ustd .> params.std_thresh] .= true
        flags[vstd .> params.std_thresh] .= true
        flags[wstd .> params.std_thresh] .= true
    end
    if params.global_validate
        flags[u .< params.u_thresh[1]] .= true
        flags[u .> params.u_thresh[2]] .= true
        flags[v .< params.v_thresh[1]] .= true
        flags[v .> params.v_thresh[2]] .= true
        flags[w .< params.w_thresh[1]] .= true
        flags[w .> params.w_thresh[2]] .= true
    end
    if params.s2n_validate
        flags[s2n .< params.s2n_thresh] .= true
    end
    return nothing
end

function replace_nans!(array, filled, max_iter, tol; kernel_size=2)
    n_dim = ndims(array)
    nan_indices = Tuple.(findall(isnan.(array)))
    nan_indices_cartesian = findall(isnan.(array))
    n_nan = length(nan_indices)
    replaced_new = zeros(Float32, n_nan)
    replaced_old = zeros(Float32, n_nan)
    for _ in 1:max_iter
        any_updates = false
        for k in 1:n_nan 
            ind = nan_indices[k]
            replaced_new[k] = 0.0
            window_bounds = [max(1, i - kernel_size):min(size(array, dim), i + kernel_size) 
                             for (dim, i) in enumerate(ind)]
            window = filled[window_bounds...]
            valid_mask = .!isnan.(window)
            if any(valid_mask)
                replaced_new[k] = sum(window[valid_mask]) / sum(valid_mask)
            else
                replaced_new[k] = NaN
            end
        end
        filled[nan_indices_cartesian] .= replaced_new
        if mean((replaced_new .- replaced_old).^2) < tol
            break
        else
            replaced_old = replaced_new
        end
    end
    return nothing
end

function replace_outliers(u, v, w, flags, params; max_iter=3, tol=1e-8, kernel_size=3)
    if params.replace_method == "localmean"
        u[flags] .= NaN
        v[flags] .= NaN 
        w[flags] .= NaN 
        u_replaced = copy(u)
        v_replaced = copy(v)
        w_replaced = copy(w)
        replace_nans!(u, u_replaced, max_iter, tol, kernel_size=kernel_size)
        replace_nans!(v, v_replaced, max_iter, tol, kernel_size=kernel_size)
        replace_nans!(w, w_replaced, max_iter, tol, kernel_size=kernel_size)
        return u_replaced, v_replaced, w_replaced 
    elseif params.replace_method == "zero"
        u[flags] .= 0.0
        v[flags] .= 0.0
        w[flags] .= 0.0
        return u, v, w
    else
        error("Invalid replace_method")
    end
end
