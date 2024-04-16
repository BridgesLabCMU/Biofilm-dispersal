using FileIO 
using NPZ
using Plots
using NaturalSort
using StatsBase
using NaNStatistics
using NaturalSort

gr()

function main()
    folder = "/mnt/h/Dispersal/WT_replicate1_processed/Displacements/"
    files = sort([f for f in readdir(folder, join=true) if occursin("trajectories", f)], lt=natural)
    _, nx, ny, nz = size(load(files[1], "trajectories"))
    trajectories = zeros(3, nx, ny, nz, length(files))
    for i in 1:length(files)
        trajectories[:,:,:,:,i] = load(files[i], "trajectories")
    end
    trajectories = trajectories[:, 1:2:end, 1:2:end, 1:2:end, 1:20]
    plt = plot(title = "Particle Trajectories", xlabel = "X", ylabel = "Y", zlabel = "Z", legend = false, xlimit=(0, nx), ylimit=(0, ny), zlimit=(0, nz))

	# Determine and plot the trajectories of moving particles
	_, nx, ny, nz, nt = size(trajectories)
    @show nx*ny*nz
	for i in 1:nx
        @show i
		for j in 1:ny
			for k in 1:nz
				x = trajectories[1, i, j, k, :]
				y = trajectories[2, i, j, k, :]
				z = trajectories[3, i, j, k, :]

				# Check if the particle has moved from its initial position
				if any(x .!= x[1]) || any(y .!= y[1]) || any(z .!= z[1])
					plot!(plt, x, y, z, label = "")  
				end
			end
		end
	end

	# Display the plot
	display(plt)
end
main()
