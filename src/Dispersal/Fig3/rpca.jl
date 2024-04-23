using LinearAlgebra
using SparseArrays
using Printf
using FileIO
using NaturalSort

function inexact_alm_rpca(D; λ=nothing, tol=1e-7, maxIter=1000)
    m, n = size(D)

    if λ === nothing
        λ = 1 / sqrt(m)
    end

    λ /= sqrt(m)

    # Initialize
    Y = copy(D)
    norm_two = norm(Y, 2)
    norm_inf = norm(Y[:], Inf) / λ
    dual_norm = max(norm_two, norm_inf)
    Y ./= dual_norm

    A_hat = zeros(m, n)
    E_hat = zeros(m, n)
    mu = 1.25 / norm_two  # Can be tuned
    mu_bar = mu * 1e7
    rho = 1.5  # Can be tuned
    d_norm = norm(D, 2)

    iter_num = 0
    converged = false
    sv = 10
    while !converged
        iter_num += 1

        temp_T = D - A_hat .+ (1 / mu) .* Y
        E_hat = max.(temp_T .- λ / mu, 0) .+ min.(temp_T .+ λ / mu, 0)

        U, S, V = svd(D - E_hat .+ (1 / mu) .* Y, full=false)
        svp = count(x -> x > 1 / mu, S)
        if svp < sv
            sv = min(svp + 1, n)
        else
            sv = min(svp + round(Int, 0.05 * n), n)
        end

        A_hat = U[:, 1:svp] * Diagonal(S[1:svp] .- 1 / mu) * V[:, 1:svp]'

        Z = D - A_hat - E_hat

        Y .+= mu .* Z
        mu = min(mu * rho, mu_bar)

        stopCriterion = norm(Z, 2) / d_norm
        if stopCriterion < tol
            converged = true
        end

        if iter_num % 10 == 0
            @printf("Iteration: %d, Rank(A): %d, ||E||_0: %d, stopCriterion: %.5f\n", 
                    iter_num, rank(A_hat), count(!iszero, E_hat), stopCriterion)
        end

        if !converged && iter_num >= maxIter
            println("Maximum iterations reached!")
            break
        end
    end

    return A_hat, E_hat, iter_num
end

function main()
    data_folders = [f for f in readdir("/mnt/h/Dispersal", join=true) if isdir(f) && occursin("processed", f)]
    for data_folder in data_folders
        data_folder = data_folder*"/Displacements/"
        @show data_folder
        files = sort([f for f in readdir(data_folder, join=true) if occursin("piv", f)], lt=natural) 
        nfiles = length(files)
        u_dummy = load(files[1], "u")
        h, width, d = size(u_dummy)
        data_matrix = Array{Float64}(undef, h*width*d*3, nfiles)
        for i in 1:nfiles
            u, v, w = load(files[i], "u", "v", "w")
            data_matrix[:,i] = vcat(u[:], v[:], w[:] .* 4) 
        end
        λ = 0.5 
        tol = 1e-7
        maxIter = 1000
        L, S, iter_num = inexact_alm_rpca(data_matrix, λ=λ, tol=tol, maxIter=maxIter)
        L = reshape(L, (h, width, d, 3, nfiles))
        @views Lu, Lv, Lw = L[:,:,:,1,:], L[:,:,:,2,:], L[:,:,:,3,:]
        save(data_folder*"rpca_result.jld2", Dict("u" => Lu, "v" => Lv, "w" => Lw))
    end
end
main()
