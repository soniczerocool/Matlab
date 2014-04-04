
result_save_dir = 'N:\User\Research\Results\BRATS-2_ICPR\';

mrf_result_save_dir =  [main_result_dir,'mrf\nodenoising\'];
mrf_denoised_result_save_dir =  [main_result_dir,'mrf\Denoised\'];

boykov_result_save_dir =  [main_result_dir,'boykov\nodenoising\'];
boykov_denoised_result_save_dir =  [main_result_dir,'boykov\Denoised\'];

bykov_ndt_result_save_dir =  [main_result_dir,'boykov_nodataterm\nodenoising\'];
bykov_ndt_denoised_result_save_dir =  [main_result_dir,'boykov_nodataterm\Denoised\'];


knn_medianfilter_result_save_dir =[main_result_dir,'knn_medianfilter\'];
knn_result_save_dir =[main_result_dir,'knn\'];

evaluation_matrix_savedir = [knn_result_save_dir,'evaluation_matrix\'];

brain_list = dir(evaluation_matrix_savedir);
brain_list = brain_list(3:end);
