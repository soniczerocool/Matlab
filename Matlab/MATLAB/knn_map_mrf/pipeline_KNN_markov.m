close all
clear all
mex -lut C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\markov_network_multilable5.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\GCoptimization.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\graph.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\LinkedBlockList.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\maxflow.cpp


%% apply graph cut
% make the graph 
% PUT T1C just for the shape of it . the actual image is not used
segmented =  markov_network_multilable5(T1C,energy,0.12);
%   segmented = u_medfilt3(segmented,5);
%segmented = approxMedfilt3(segmented,5,5,5)
%% assign labels 
[h,w,d] = size(T1C);
segmented_volume = label_to_matrix(segmented,space,h,w,d);

time = toc(t);
%% compute statistics 
segmented_brain = zeros(height,width,depth);
segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;

%[dice,jaccard,precision,recall] = compute_eval_metrics(segmented,truth);


%% save results
% 
%% save in nifti format

segmented_nii = make_nii(segmented); 
save_nii(segmented_nii, [result_path_brain,'\',type,'_',name,'_segmented_tumor.nii.gz']);

 %% save in mha format
 mha_result_path = 'markov_segment_beta_0.12222.mha';
segmented_brain=uint16(segmented_brain);
 writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset);

