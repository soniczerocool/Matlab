close all
 clear all
mex -lut C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\markov_network_multilable_general3Dgraph.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\GCoptimization.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\graph.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\LinkedBlockList.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\maxflow.cpp
%%
global T1C T2 FLAIR
%% define paths and parameters
Nb_Selected = 500;
k=5;
Nb_healthy = 500;
downsamplerate = 20;
 results_path_root = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\CHALLENGE2013_KNN_mha_test\';
 type = 'HG';
 name='0001';
 Brain_path = ['N:\User\Research\Brain_data\BRATS-2\Image_Data\',type,'\',name];
 if exist(results_path_root,'dir')==0
    mkdir(results_path_root)
end
result_path_brain = [results_path_root,type,'\', name];
if exist(result_path_brain,'dir')==0
    mkdir(result_path_brain)
end
mask_name_path = 'N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\HG\0001';


%% load modalities
[T1C, T2, FLAIR, truth, info] = load_modalities(Brain_path);
[height width depth] = size(T1C);
T1C = double(T1C);
T2 = double(T2);
FLAIR = double(FLAIR);

%% creating the mask
%MASK = make_mask(T1C, T2, FLAIR);
% or load the mask from a path (future work)
 MASK = load_mask(mask_name_path);

 [ymin,ymax,xmin,xmax,zmin,zmax] = make_regtangle(T1C);
T1C = T1C(xmin:xmax , ymin:ymax , zmin:zmax);
T2 = T2(xmin:xmax , ymin:ymax , zmin:zmax);
FLAIR = FLAIR(xmin:xmax , ymin:ymax , zmin:zmax);

MASK = MASK(xmin:xmax , ymin:ymax , zmin:zmax);
truth = truth(xmin:xmax , ymin:ymax , zmin:zmax);
truth_backup = truth;
%% make space and selected points space
[space,mask_idx,truth_idx] = make_space(T1C, FLAIR, T2,truth, MASK);

[selected_space, mask_idx] = make_selected_space(space, mask_idx, Nb_healthy,Nb_Selected );

%% KNN search

%importance_map = KNN_search(space, selected_space, k);
%[IDX,D] = KNN_search(selected_space,space, k)
[IDX,D] = knnsearch(selected_space',space','K', k ,'NSMethod','kdtree' ,'Distance','euclidean');

%% assign labels
matrix = make_dataterm_matrix(IDX,mask_idx,space);
%normalize it
matrix = (matrix/k);
% convert to energy
%energy = -log(matrix);
energy = 1 - matrix;
t=tic;

%% apply graph cut
% make the graph 
% PUT T1C just for the shape of it . the actual image is not used
segmented =  markov_network_multilable_general3Dgraph(T1C,energy,0.2);
%   segmented = u_medfilt3(segmented,5);
%segmented = approxMedfilt3(segmented,5,5,5)
%% assign labels 
[h,w,d] = size(T1C);
segmented_volume = label_to_matrix(segmented,space,h,w,d);

%time = toc(t);
%% compute statistics 
segmented_brain = zeros(height,width,depth);
segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
truth = truth_backup;
[dice,jaccard,precision,recall] = compute_eval_metrics(segmented_brain,truth);


%% save results
% 
% %save nii files
% mask_nii = make_nii(MASK);
% save_nii(nii, [result_path_brain,'\','6dimportancemap.nii.gz']);
% save_nii(mask_nii, [result_path_brain,'\','MASK.nii.gz']);
 %segmented_nii = make_nii(segmented); 
 %save_nii(segmented_nii [result_path_brain,'\',type,'_',name,'_segmented_tumor.nii.gz']);
% 
% %save statistics
% save([result_path_brain,'\',type,'_',name,'_stat_results.mat'],'precision','recall');
% 
% %save meta files(.mha)
%mha_result_name = strrep(flair_name , 'XX.O.MR_Flair',['Seg_',type,'_',name]);
mha_result_path = 'markov_segment_beta_0.122248.mha';
segmented_brain=uint16(segmented_brain);
 writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset);
% MASK = uint16(MASK);
% mha_MASK_path = [result_path_brain,'\','MASK_',type,'_',name,'.mha'];
% writemetaimagefile(mha_MASK_path, MASK, infot1c.PixelDimensions,infot1c.Offset);
% 
