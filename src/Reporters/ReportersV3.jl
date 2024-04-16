using Images: imfilter, mapwindow, adjust_histogram!, LinearStretching, Gray, N0f16, N0f8, Kernel, RGB 
using ImageView
using HistogramThresholding: find_threshold, Otsu
using StatsBase: mean, median, quantile
using TiffImages
using FileIO
using NaturalSort: sort, natural
using FLoops

function calculate_background!(stack, background, ntimepoints)
    for t in 1:ntimepoints
        @views background[:,:,t] = imfilter(stack[:,:,t], Kernel.gaussian(50))
    end
    return nothing
end

function compute_mask!(stack, background, masks, ntimepoints, constitutive)
    imgs_noback = similar(stack)
    @inbounds for t in 1:ntimepoints
        @views img = imfilter(stack[:,:,t] - background[:,:,t], Kernel.gaussian(2)) 
        imgs_noback[:,:,t] = img
	end
    thresh = find_threshold(imgs_noback, Otsu()) - 0.00037
    @inbounds for t in 1:ntimepoints
        @views mask = imgs_noback[:,:,t] .> thresh
		masks[:,:,t] = mask
    end
end

function output_images!(stack, masks, overlay, output_dir, replicate)
    imshow(masks)
    stack .*= masks
    flat_stack = vec(stack)
    img_min = quantile(flat_stack, 0.0035)
    img_max = quantile(flat_stack, 0.9965)
    adjust_histogram!(stack, LinearStretching(src_minval=img_min, src_maxval=img_max, 
                                              dst_minval=0, dst_maxval=1))
    TiffImages.save("$output_dir/constitutive_"*replicate*".tif", stack)
    @inbounds for i in CartesianIndices(stack)
        gray_val = RGB{N0f8}(stack[i], stack[i], stack[i])
        overlay[i] = masks[i] ? RGB{N0f8}(1, 0, 0) : gray_val
    end
    TiffImages.save("$output_dir/constitutive_mask_"*replicate*".tif", overlay)
end

function main()
    dir = "/mnt/f/Sandhya_Imaging/Time_Lapses/Reporters/vpsL/" # Folder containing the images for strain X 
    files = sort([f for f in readdir(dir, join=true) if occursin(".ome.tif", f)], lt=natural)
    sig = 2 
    constitutive = 2 # index for constitutive channel
    reporter = 1 # index for reporter channel
    ntimepoints = 17
    if isdir("$dir/results_images")
        rm("$dir/results_images"; recursive = true)
    end
    if isdir("$dir/results_data")
        rm("$dir/results_data"; recursive = true)
    end
    output_img_dir = "$dir/results_images"
    output_data_dir = "$dir/results_data"
    mkdir(output_img_dir)
    mkdir(output_data_dir)
    output_file = "$output_data_dir/data.jld2"
    data_matrix = Array{Float64, 2}(undef, ntimepoints, length(files))
    for (j, file) in enumerate(files)
        images = TiffImages.load(file)
        height, width, _ = size(images)
        metadata = first(ifds(images))[TiffImages.IMAGEDESCRIPTION].data
        Zloc = findfirst("SizeZ=", metadata)
        Cloc = findfirst("SizeC=", metadata)
        Tloc = findfirst("SizeT=", metadata)
        nslices = tryparse(Int, metadata[(Zloc[4]+4):(Zloc[4]+5)])
        channels = tryparse(Int, string(metadata[Cloc[4]+4]))
        nframes = tryparse(Int, metadata[(Tloc[4]+4):(Tloc[4]+5)])
        @views images = sum(Float64.(reshape(images, (height, width, nslices, channels, nframes))), dims=3)[:,:,1,:,:]
        TiffImages.save("$output_img_dir/test.tif", Gray{N0f16}.(images[:,:,constitutive,:])) 
        @views output_stack = images[:,:,constitutive,:]
        masks = zeros(Bool, size(output_stack))
        background = similar(output_stack)
        calculate_background!(output_stack, background, nframes)
        compute_mask!(output_stack, background, masks, 
                      nframes, constitutive)
        @views reporter_stack = images[:,:,reporter,:]
        let output_stack = output_stack 
            @floop for t in 1:nframes
                @inbounds signal = @views mean((reporter_stack[:,:,t] ./ output_stack[:,:,t]) .* masks[:,:,t])
                @inbounds data_matrix[t, j] = signal 
            end
        end
        output_stack = Gray{N0f8}.(output_stack)
        overlay = zeros(RGB{N0f8}, size(output_stack)...)
        output_images!(output_stack, masks, overlay, output_img_dir, string(j))
    end
    data_matrix .= ifelse.(isnan.(data_matrix), 0, data_matrix)
    FileIO.save(output_file, Dict("data" => data_matrix))
end
main()
