using FileIO 
using NPZ
using NaturalSort
using StatsBase
using NaNStatistics
using NaturalSort
using Revise 
using GLMakie

function main()
    folder = "/mnt/h/Dispersal/WT_replicate1_processed/Displacements/"
    files = sort([f for f in readdir(folder, join=true) if occursin("particle_state", f)], lt=natural)
    _, nx, ny, nz = size(load(files[1], "trajectories"))
    trajectories = zeros(3, nx, ny, nz, length(100:450))
    for i in 1:length(100:450)
        trajectories[:,:,:,:,i] = load(files[i], "trajectories")
    end
    trajectories = trajectories[:, 1:3:end, 1:3:end, 1:3:end, :]
    plt = GLMakie.Figure(fontsize = 30)
    Axis3(plt[1, 1])

	# Determine and plot the trajectories of moving particles
	_, nx, ny, nz, nt = size(trajectories)
    t = 1:nt
    @show nx*ny*nz
	for i in 1:nx
        @show i
		for j in 1:ny
			for k in 1:nz
				x = trajectories[1, i, j, k, :]
				y = trajectories[2, i, j, k, :]
				z = trajectories[3, i, j, k, :]

				# Check if the particle has moved from its initial position
                dist = sqrt((x[end] - x[1])^2 + (y[end] - y[1])^2 + (z[end] - z[1])^2)
                if dist > 20 && dist < 100 && !any(<(0), z) 
                    ps = Point3f[[xi, yi, zi] for (xi, yi, zi) in zip(x, y, z)]
                    @show ps
					lines!(ps, color=t, linewidth=2, colormap=:inferno)  
				end
			end
		end
	end

	# Display the plot
	plt
end
main()
