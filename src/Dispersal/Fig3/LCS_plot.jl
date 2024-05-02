using GLMakie
using FileIO
using NaturalSort
using LinearAlgebra
using Images 
using PythonCall
using StatsBase
using ImageMorphology

filters = pyimport("skimage.filters")


function set_boundary_voxels_to_zero(image, distance)
    nx, ny, nz = size(image)
    mask = trues(size(image))
    mask[distance+1:nx-distance, distance+1:ny-distance, 1:nz-distance] .= false
    dist_transform = distance_transform(feature_transform(mask))
    image[dist_transform .<= distance] .= false
    return image
end

function extract_lcs(ftle_field, thresh)
    sato_filtered = pyconvert(Array, filters.sato(Py(ftle_field).to_numpy(), [3], black_ridges=false))
    LCS = sato_filtered .> thresh 
	LCS = set_boundary_voxels_to_zero(LCS, 20)
    return LCS
end

function main()
    # Load data
    folder = "/mnt/h/Dispersal/WT_replicate1_processed/Displacements/"
    images_folder = "/mnt/h/Dispersal/WT_replicate1_processed/"
    mask_files = sort([f for f in readdir(images_folder, join=true) if occursin("mask_isotropic", f)], lt=natural)
    mask_t1 = load(mask_files[10])
    files = sort([f for f in readdir(folder, join=true) if occursin("FTLE", f)], lt=natural)
    piv_files = sort([f for f in readdir(folder, join=true) if occursin("piv", f)], lt=natural)
    dummy_file = load(files[1], "forward_LCS")
    height, width, depth = size(dummy_file)
    forward_LCS = zeros(Float32, height, width, depth, length(files))
    backward_LCS = zeros(Float32, height, width, depth, length(files))
    for i in 1:length(files)-1
        forward_LCS[:,:,:,i], backward_LCS[:,:,:,i] = load(files[i], "forward_LCS", "backward_LCS")
    end
    backward_plot = sum(backward_LCS[:,:,:,10:45], dims=4)[:,:,:,1]
    backward_LCS = extract_lcs(backward_plot, 0.2)
    fig = Figure(fontsize=30)
    Axis3(fig[1, 1])
    ps = Point3f[[xi, yi, zi] for xi in 1:height, yi in 1:width, zi in 1:depth if backward_LCS[xi, yi, zi] == 1]
    cs = [zi for xi in 1:height, yi in 1:width, zi in 1:depth if backward_LCS[xi, yi, zi] == 1]
    scatter!(ps, color=cs, colormap=:magma, markersize = 20)
    fig
end
main()
