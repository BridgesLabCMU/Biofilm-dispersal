using FileIO 
using NPZ
using Plots
using NaturalSort
using StatsBase
using NaNStatistics

function main()
    folder = "/mnt/h/Dispersal/WT_replicate1_processed/Displacements/"
    file = [f for f in readdir(folder, join=true) if occursin("trajectories", f)][1]

    trajectories = load(file, "trajectories") # (3, nx, ny, nz, nt) -- (nx, ny, nz) indexes each particle 
    _, nx, ny, nz, nt = size(trajectories)
	displacements = sum(abs.(trajectories[:, :, :, :, 1:nt] .- trajectories[:, :, :, :, 1]), dims=1)

	# Find indices where displacement is non-zero
	moving_particles_indices = findall(x -> x > 0, displacements)

	# Prepare data for plotting
	xs, ys, zs, times = [], [], [], []
	for ind in moving_particles_indices
		# Extracting the trajectory of the moving particle
		traj = trajectories[:, ind[2], ind[3], ind[4], :]
		for t in 1:nt
			push!(xs, traj[1, t])
			push!(ys, traj[2, t])
			push!(zs, traj[3, t])
			push!(times, t)
		end
	end

    @show size(xs)

	# Plotting
	pgfplotsx() # Interactive 3D plot
	scatter(xs, ys, zs, zcolor=times, palette=:viridis, legend=false, xlabel="X", ylabel="Y", zlabel="Z", title="Particle Trajectories")

    #save("$folder/quiver_test_inverted.png", fig)
end
main()
