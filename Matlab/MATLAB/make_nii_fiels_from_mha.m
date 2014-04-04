% make niftifiles from mha

seg_path =   'N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\LG\0004\segmented_LG_0004.mha'
truth_path = 'N:\User\Research\Brain_data\BRATS-2\Image_Data\LG\0004\VSD.Brain_2more.XX.XX.OT\VSD.Brain_2more.XX.XX.OT.6604.mha'

output_path_seg = 'L:\Proj_intctiv_brain_Seg\nii_mha_files_comparision\LG_0004_seg.nii.gz'
output_path_truth = 'L:\Proj_intctiv_brain_Seg\nii_mha_files_comparision\LG_0004_truth.nii.gz'

seg =  mha_read_volume(seg_path);
truth = mha_read_volume(truth_path);

%IMAGE = load_nii(output_path, [], 1);
%data =IMAGE.img;

seg_nii = make_nii(seg);
save_nii(seg_nii, output_path_seg);

truth_nii = make_nii(truth);
save_nii(truth_nii, output_path_truth);


[dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(seg,truth)