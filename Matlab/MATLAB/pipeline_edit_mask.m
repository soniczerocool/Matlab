clear all
close all
global T1C T2 FLAIR

%% define paths and parameters
Nb_Selected = 500;
k=21;
Nb_healthy = 500;
downsamplerate = 20;
 results_path_root = 'N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\';
%results_path_root = '/Users/uoft/Dropbox/PhD/BratsAnalysis/CHALLENGE2013_KNN_mha/';
 type = 'HG';
 name='0004';
 Brain_path = ['N:\User\Research\Brain_data\BRATS-2\Image_Data\',type,'\',name];
% Brain_path = ['/Volumes/FreeAgent Drive/Research UdeSh/BRATS-2/Image_Data/',type,'/',name];
 mask_path = ['N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\',type,'\',name]
 if exist(results_path_root,'dir')==0
    mkdir(results_path_root)
end
result_path_brain = [results_path_root,type,'\', name];
if exist(result_path_brain,'dir')==0
    mkdir(result_path_brain)
end

%% load modalities
[T1C, T2, FLAIR, truth, info] = load_modalities(Brain_path);
[height width depth] = size(T1C);
T1C = double(T1C);
T2 = double(T2);
FLAIR = double(FLAIR);

%% creating the mask
%MASK = make_mask(T1C, T2, FLAIR);
% or load the mask from a path (future work)
 MASK = load_mask(mask_path);
 MASK = double(MASK);
MASK = edit_mask(T1C, T2, FLAIR,MASK);





% %% assign labels of the mask form truth
% % putting the healthy labels to 10
% h_index = find(truth==0);
% truth(h_index) = 10;
% 
% 
% Background_index = find(MASK > 0);
% MASK(Background_index) = truth(Background_index);
% nicroses_idx = find(MASK ==1);
% MASK(nicroses_idx) = 3; 

%% make space and selected points space
%[space,mask_idx] = make_space(T1C, FLAIR, T2,truth, MASK);

%[selected_space, mask_idx] = make_selected_space(space, mask_idx, Nb_healthy,Nb_Selected );

%% KNN search

%importance_map = KNN_search(space, selected_space, k);
%[IDX,D] = KNN_search(selected_space,space, k)
%[IDX,D] = knnsearch(selected_space',space','K', k ,'NSMethod','kdtree' ,'Distance','euclidean');

%% assign labels
%segmented = assign_labels(IDX,D,mask_idx, space, height, width, depth);

%t=tic;

%% apply median filter
%segmented = medfilt3(segmented, [5,5,5]);
%   segmented = u_medfilt3(segmented,5);
%time = toc(t);
%% compute statistics 
 [dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(MASK,truth);


%% save results
% save the mha file for the mask
MASK = uint16(MASK);
mha_MASK_path = [result_path_brain,'\','MASK_',type,'_',name,'.mha'];
writemetaimagefile(mha_MASK_path, MASK, info.PixelDimensions,info.Offset);


%save nii files
%mask_nii = make_nii(MASK);
%save_nii(nii, [result_path_brain,'\','6dimportancemap.nii.gz']);
%save_nii(mask_nii, [result_path_brain,'\','MASK.nii.gz']);
%segmented_nii = make_nii(segmented); 
%save_nii(segmented_nii, [result_path_brain,'\',type,'_',name,'_segmented_tumor.nii.gz']);

%%save statistics
%save([result_path_brain,'\',type,'_',name,'_stat_results.mat'],'precision','recall');

%%save meta files(.mha)
%mha_result_name = strrep(flair_name , 'XX.O.MR_Flair',['Seg_',type,'_',name]);
%mha_result_path = [result_path_brain,'\',mha_result_name];
%segmented=uint16(segmented);
%writemetaimagefile(mha_result_path, segmented, info.PixelDimensions,info.Offset);
