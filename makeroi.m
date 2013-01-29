function mask = makeroi(struct_mask_id_arr, seg_file, outfile)

    segmentation_file = MRIread(seg_file);
    segmentation_volume = segmentation_file.vol;
    mask_volume = zeros(size(segmentation_volume));
    for i = 1:length(struct_mask_id_arr)
        mask_volume = mask_volume + (segmentation_volume ==struct_mask_id_arr(i));
    end
    mask = segmentation_file;
    mask.vol = mask_volume;
    MRIwrite(mask,outfile);

end

