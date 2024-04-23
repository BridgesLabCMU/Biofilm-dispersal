using GLMakie
using FileIO

function main()
    # Load data
    folder = "/mnt/h/Dispersal/WT_replicate1_processed/Displacements/"
    files = [f for f in readdir(folder, join=true) if occursin("FTLE", f)]
    dummy_file = load(files[1], "forward_LCS")
    height, width, depth = size(dummy_file)
    forward_LCS = zeros(Float32, height, width, depth, length(files))
    backward_LCS = zeros(Float32, height, width, depth, length(files))
    for i in 1:length(files)
        forward_LCS[:,:,:,i], backward_LCS[:,:,:,i] = load(files[i], "forward_LCS", "backward_LCS")
    end
    forward_plot = sum(forward_LCS[:,:,:,15:55], dims=4)[:,:,:,1]
    backward_plot = sum(backward_LCS[:,:,:,15:55], dims=4)[:,:,:,1]
    fig = Figure(size=(1000, 450), fontsize = 30)
    colormap = to_colormap(:plasma)
    colormap[1] = RGBAf(0,0,0,0)
    volume(fig[1, 1], forward_plot, algorithm = :iso, colormap=colormap, axis=(type=Axis3, title = "Forward FTLE"))
    volume(fig[1, 2], backward_plot, algorithm = :iso, colormap=colormap, axis=(type=Axis3, title = "Backward FTLE"))
    fig
end
main()
