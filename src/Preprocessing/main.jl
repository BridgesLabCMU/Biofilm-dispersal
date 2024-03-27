include("CleanTimelapse.jl")
using .CleanTimelapse

cell_threshold = 3.0e-5 # Depends on imaging setup!!
file_directory = "/mnt/h/Dispersal/test/"
zstack_thresh = 7.5 
clean_timelapse(cell_threshold, file_directory, zstack_thresh)
