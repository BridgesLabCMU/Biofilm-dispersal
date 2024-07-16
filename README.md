Code for the manuscript "Biofilm dispersal dynamics revealed using far-red fluorogenic probes"

The code is structured by relevant figure. In the source folder, the code related to all figures pertaining to dispersal is contained in the "Dispersal" folder, while the code used to analyze the results in Fig. S2 and Fig. S3 is contained in distinct folders. Code for preprocessing the high-resolution timelapses is contained in the "Preprocessing" folder. Code for computing intra-biofilm displacements is contained in the "Image-cross-correlation" folder.

To analyze high-resolution timelapses, the timelapses were first preprocessed using the code contained in the "Preprocessing" folder. This can be done by simply editing and running the "main.jl" file. Some lines related to file name strucutre may need to be edited in the "CleanTimelapse.jl" file. Running the code will generate n registered z-stacks and n registered z-stacks with planktonic cells omitted, where n is the number of timepoints, in a new folder. Moving
forward, the registered z-stacks with planktonic cells omitted are used. Each step in the subsequent analysis can then be performed by running a relevant script in either the "Dispersal" folder or the "Image-cross-correlation" folder. Displacements, for example, can be computed by editing paths and running the "main.jl" file in the "Image-cross-correlation" folder.

If you find bugs, please contact jojo@cmu.edu.
