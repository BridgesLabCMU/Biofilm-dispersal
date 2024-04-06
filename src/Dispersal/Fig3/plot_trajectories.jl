using FileIO 
using NPZ
using CairoMakie
using NaturalSort
using StatsBase
using NaNStatistics

function main()
    folder = "/mnt/h/Dispersal/WT_replicate1_processed/Displacements/"
    file = [f for f in readdir(folder, join=true) if occursin("trajectories", f)][1]

    trajectories = load(file, "trajectories")

    save("$folder/quiver_test_inverted.png", fig)
end
main()
