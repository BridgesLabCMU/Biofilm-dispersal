function fillWindow!(window, image, i, j, k, overlap)
    i1 = i-overlap[1]
    i2 = i+overlap[1]-1
    j1 = j-overlap[2]
    j2 = j+overlap[2]-1
    k1 = k-overlap[3]
    k2 = k+overlap[3]-1
    for x in i1:i2 
        for y in j1:j2 
            for z in k1:k2 
                window[y-j1+1, x-i1+1, z-k1+1] = complex(image[y, x, z]) 
            end
        end
    end
    return nothing
end

function cross_correlate!(corr, window2, window1, plan, iplan)
    plan * window2
    plan * window1
    @inbounds @simd for i in 1:length(window2)
        window2[i] = window2[i] * conj(window1[i])
    end
    iplan * window2
    @inbounds @simd for i in 1:length(window2)
        corr[i] = real(window2[i]) 
    end
    return nothing
end

Memoize.@memoize function firstPeak(corr)
	val, idx = findmax(corr)
    return val, Tuple(idx)
end

function gaussian_subpixel(neighbors, maxval)
    up = log(neighbors[13])
    down = log(neighbors[15])
    left = log(neighbors[11])
    right = log(neighbors[17])
    front = log(neighbors[5])
    back = log(neighbors[23])
    mx = log(maxval)
	return [(up-down)/(2*up-4*mx+2*down),
             (left-right)/(2* left-4*mx+2*right),
             (front-back)/(2*front-4*mx+2*back)];
end

neighboridx(j) = [(j[1]+y, j[2]+x, j[3]+z) for y in -1:1, x in -1:1, z in -1:1]

function corr_to_disp!(corr, U, V, W, grid_idx)
    corr_center = div.(size(corr), 2)
    val, peak = firstPeak(corr)
    if any(i == 1 for i in peak) || any(peak .== size(corr)) 
        disp = peak .- corr_center .- 1
        U[grid_idx] = disp[2]
        V[grid_idx] = disp[1] 
        W[grid_idx] = disp[3] 
    else
        neighbors = [corr[j...] for j in neighboridx(peak)]
        minVal = minimum(neighbors) 
        for j in 1:length(neighbors)
            neighbors[j] = 1 + neighbors[j] - minVal
        end
        refinement = gaussian_subpixel(neighbors, 1+val-minVal)
        refinement = [isnan(x) ? 0.0 : x for x in refinement] 
        disp = peak .+ refinement .- corr_center .- 1
        U[grid_idx] = disp[2]
        V[grid_idx] = disp[1]
        W[grid_idx] = disp[3]
    end
    return nothing  
end

function compute_snr!(corr, SNR, grid_idx, method) 
    max2 = 0.0 
    h, w, d = size(corr,1), size(corr,2), size(corr,3)
    max1, idx1 = firstPeak(corr) 
    empty!(memoize_cache(firstPeak))
    if idx1[1] == 1 || idx1[1] == h || idx1[2] == 1 || 
                        idx1[2] == w || idx1[3] == 1 || idx1[3] == d
        SNR[grid_idx] = 0.0
    else
        if method == "peak2mean"
            max2 = mean(corr)
        end
        if max2 < 0 && max1 > 0 
            max1 = max1 - max2 + 1
            max2 = 1
        end
        SNR[grid_idx] = max2 !== 0.0 ? max1 / max2 : 0.0
    end
    return nothing
end

function compute_displacements!(image1, image2, window_size, x, y, z, 
                                overlap, s2n_method, U, V, W, SNR)
    # Pre-allocate
    window1 = zeros(Complex{Float32}, window_size[1], 
                    window_size[2], window_size[3])
    window2 = zeros(Complex{Float32}, window_size[1], 
                    window_size[2], window_size[3])
    corr = zeros(Float32, window_size[1], window_size[2], window_size[3])
    plan = FFTW.plan_fft!(window2)  
    iplan = FFTW.plan_ifft!(window2)
    @inbounds for k in z 
        @inbounds for i in x 
            @inbounds for j in y 
                grid_idx = CartesianIndex([j ÷ overlap[1], i ÷ overlap[2], k ÷ overlap[3]]...)
                window1 .= 0
                window2 .= 0
                fillWindow!(window1, image1, i, j, k, overlap)
                fillWindow!(window2, image2, i, j, k, overlap)
                cross_correlate!(corr, window2, window1, plan, iplan)
                corr = FFTW.fftshift(corr)
                corr_to_disp!(corr, U, V, W, grid_idx)
                compute_snr!(corr, SNR, grid_idx, s2n_method)
            end
        end
    end
    return nothing 
end
