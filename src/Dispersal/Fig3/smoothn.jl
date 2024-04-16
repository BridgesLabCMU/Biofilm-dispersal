using FFTW
using LinearAlgebra
using Statistics
using Distributions
using Optim
using FileIO

function smoothn(y; nS0=10, axis=nothing, smoothOrder=2.0, sd=nothing, verbose=false, s0=nothing, z0=nothing, isrobust=false, W=nothing, s=nothing, MaxIter=100, TolZ=1e-3, weightstr="bisquare")
    if sd !== nothing
        W = zeros(size(sd))
        mask = sd .> 0
        W[mask] .= 1.0 ./ sd[mask] .^ 2
    end
    if W !== nothing
        W ./= maximum(W) 
    end

    sizy = size(y)
    if axis === nothing 
        axis = 1:ndims(y)
    end
    noe = length(y)

    if noe < 2
        return y, s, 0, W
    end
    if W === nothing
        W = ones(sizy) 
    end
    weightstr = lowercase(weightstr)
    IsFinite = isfinite.(y)
    nof = sum(IsFinite)
    W .= W .* IsFinite
    isweighted = any(W .!= 1)
    isauto = s === nothing

    Lambda = zeros(sizy)
    for i in axis
        siz0 = ones(Int, ndims(y))
        siz0[i] = size(y, i)
        Lambda .+= reshape(cos.((collect(1:siz0[i]) .- 1.0) .* pi / siz0[i]), Tuple(siz0))
    end
    Lambda = -2.0 .* (length(axis) .- Lambda)

    if !isauto
        Gamma = 1.0 ./ (1 .+ (s .* abs.(Lambda)) .^ smoothOrder)
    end

    N = sum(size(y) .!= 1)
    hMin = 1e-6
    hMax = 0.99
    sMinBnd = sqrt((((1 + sqrt(1 + 8 * hMax ^ (2.0 / N))) / 4 / hMax ^ (2.0 / N)) ^ 2 - 1) / 16.0)
    sMaxBnd = sqrt((((1 + sqrt(1 + 8 * hMin ^ (2.0 / N))) / 4 / hMin ^ (2.0 / N)) ^ 2 - 1) / 16.0)

    Wtot = W
    z = ifelse(isweighted, z0 !== nothing ? z0 : copy(y), zeros(sizy))
    z0 = copy(z)
    y[.!IsFinite] .= 0
    tol = 1.0
    RobustIterativeProcess = true
    RobustStep = 1
    nit = 0
    RF = 1 + 0.75 * isweighted

    exitflag = false
    while RobustIterativeProcess
        aow = sum(Wtot) / noe
        while tol > TolZ && nit < MaxIter
            if verbose
                println("tol: ", tol, " nit: ", nit)
            end
            nit += 1
            DCTy = dctnd(Wtot .* (y .- z) .+ z)
            if isauto && !isinteger(log2(nit))
                if s0 === nothing
                    ss = collect(LinRange(log10(sMinBnd), log10(sMaxBnd), nS0))
                    g = [gcv(p, Lambda, aow, DCTy, IsFinite, Wtot, y, nof, noe, smoothOrder) for p in ss]
                    xpost = ss[argmin(g)]
                else
                    xpost = log10(s0)
                end
                res = optimize(p -> gcv(10^p[1], Lambda, aow, DCTy, IsFinite, Wtot, y, nof, noe, smoothOrder), [log10(sMinBnd)], [log10(sMaxBnd)], [xpost], Fminbox())
                xpost = Optim.minimizer(res)[1]
                s = 10^xpost
                s0 = xpost
            end

            Gamma = 1.0 ./ (1 .+ (s .* abs.(Lambda)) .^ smoothOrder)
            z = RF * dctnd(Gamma .* DCTy) .+ (1 - RF) .* z
            tol = isweighted * norm(z0 - z) / norm(z)
            z0 = copy(z)
        end
        exitflag = nit < MaxIter

        if isrobust
            h = sqrt(1 + 16.0 * s)
            h = sqrt(1 + h) / sqrt(2) / h
            h = h ^ N
            Wtot = W .* robustweights(y - z, IsFinite, h, weightstr)
            isweighted = true
            tol = 1
            nit = 0
            RobustStep += 1
            RobustIterativeProcess = RobustStep < 3
        else
            RobustIterativeProcess = false
        end
    end

    if isauto
        if abs(log10(s) - log10(sMinBnd)) < 0.1
            println("Warning: the lower bound for s has been reached.")
        elseif abs(log10(s) - log10(sMaxBnd)) < 0.1
            println("Warning: the upper bound for s has been reached.")
        end
    end

    return z, s, exitflag, Wtot
end

function robustweights(r, I, h, wstr)
    MAD = median(abs.(r[I] .- median(r[I])))
    u = abs.(r ./ (1.4826 * MAD) ./ sqrt(1 .- h))
    c = 0.0
    W = zeros(length(r))

    if wstr == "cauchy"
        c = 2.385
        W .= 1.0 ./ (1 .+ (u ./ c) .^ 2)
    elseif wstr == "talworth"
        c = 2.795
        W .= u .< c
    else  
        c = 4.685
        W .= (1 .- (u ./ c) .^ 2) .^ 2 .* (u ./ c .< 1)
    end

    replace!(W, NaN=>0)
    return W
end

function gcv(p, Lambda, aow, DCTy, IsFinite, Wtot, y, nof, noe, smoothOrder)
    s = 10^p
    Gamma = 1.0 ./ (1 .+ (s .* abs.(Lambda)) .^ smoothOrder)
    if aow > 0.9
        RSS = norm(DCTy .* (Gamma .- 1.0))^2
    else
        yhat = dctnd(Gamma .* DCTy)
        RSS = norm(sqrt.(Wtot[IsFinite]) .* (y[IsFinite] .- yhat[IsFinite]))^2
    end
    TrH = sum(Gamma)
    GCVscore = RSS / nof / (1.0 - TrH / noe)^2
    return GCVscore
end

function dctnd(data)
    if ndims(data) == 1
        return dct(data)
    elseif ndims(data) == 2
        return dct(dct(data, 1), 2)
    elseif ndims(data) == 3
        return dct(dct(dct(data, 1), 2), 3)
    end
end

function main()
    vector_folder = "/mnt/h/Dispersal/WT_replicate1_processed/Displacements/"
    vector_files = [f for f in readdir(vector_folder, join=true) if occursin("piv", f)]
    for (i, file) in enumerate(vector_files) 
        u, v, w = load(file, "u", "v", "w")
        u, d1, d2, d3 = smoothn(u, s=0.01)
        v, d1, d2, d3 = smoothn(v, s=0.01)
        w, d1, d2, d3 = smoothn(w, s=0.01)
        save(vector_folder*"smoothed_$(i).jld2", Dict("u" => u, "v" => v, "w" => w))
    end
end

main()
