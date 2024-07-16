include("CleanTimelapse.jl")
using .CleanTimelapse

cell_threshold = 3.0e-5 # Depends on imaging setup!!
file_directory = "/mnt/f/Sandhya_Imaging/Software/"
zstack_thresh = 7.5 
CleanTimelapse.clean_timelapse(cell_threshold, file_directory, zstack_thresh)
