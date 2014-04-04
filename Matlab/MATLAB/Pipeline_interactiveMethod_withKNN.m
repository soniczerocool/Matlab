clear all
close all
global T1C T2 FLAIR
%% define paths and parameters
Nb_Selected = 500;
k=5;
Nb_healthy = 500;
downsamplerate = 20;
 results_path_root = 'N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\HG\0001\2';
 type = 'HG';
 name='0001';
  Brain_path = ['N:\User\Research\Brain_data\BRATS-2\Image_Data\',type,'\',name];
  %Brain_path = 'C:\Users\havm2701\Dropbox\TumorFortinMohammad\';
 if exist(results_path_root,'dir')==0
    mkdir(results_path_root)
end
result_path_brain = [results_path_root,type,'\', name];
if exist(result_path_brain,'dir')==0
    mkdir(result_path_brain)
end

%% load modalities
[T1C, T2, FLAIR, truth, info,flair_name] = load_modalities(Brain_path);
%[T1C, T2,FLAIR] = load_modalities_fortin(Brain_path);
[height width depth] = size(T1C);
T1C = double(T1C);
T2 = double(T2);
%FLAIR = zeros(size(T1C));
FLAIR = double(FLAIR);
%truth = 0;


%% creating the mask
%MASK = make_mask(T1C, T2, FLAIR);
% or load the mask from a path (future work)
mask_name_path = 'N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\HG\0001';
 MASK = load_mask(mask_name_path);

%% make space and selected points space
[space,mask_idx,truth_idx]  = make_space(T1C, FLAIR, T2,truth, MASK);



background = find(sum(space(1:3,:))< 0.0001);
 space(:,background) = [];
  mask_idx(background) = [];

[selected_space, mask_idx] = make_selected_space(space, mask_idx, Nb_healthy,Nb_Selected );

%make_txt_file(space,selected_space,truth_idx,mask_idx,type, name);
%% remove spatial information

%selected_space1 = selected_space(4:6,:);
%space1 = space(4:6,:);
%selected_space = selected_space(1:3,:);
%space = space(1:3,:);
%[IDX,D] = KNN_search(selected_space,space, k);
%% KNN search

%[IDX,D] = knnsearch(selected_space',space','K', k ,'NSMethod','kdtree' ,'Distance','euclidean');

[ind,knndist] = kdtreeKNN(selected_space,mask_idx,space,k);
IDX = ind';
D = knndist';
%% putback spatial information

%selected_space = [selected_space ; selected_space1];
%space = [ space;space1];
%% assign labels
segmented = assign_labels(IDX,D,mask_idx, space, height, width, depth);

t=tic;

%% apply median filter
segmented = medfilt3(segmented, [5,5,5]);
%   segmented = u_medfilt3(segmented,5);
%segmented = approxMedfilt3(segmented,5,5,5)
time = toc(t);
%% compute statistics 
%[dice,jaccard,precision,recall] = compute_eval_metrics(segmented,truth);

[dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented, truth)
%% save results

%save nii files
 %mask_nii = make_nii(MASK,[1 1 1],[8.3895 21.5912  25.1251] , 4);
 %save_nii(nii, [result_path_brain,'\','6dimportancemap.nii.gz']);
 %save_nii(mask_nii, ['C:\Users\havm2701\Dropbox\Tumor_fortin\','MASK3.nii.gz']);
 %segmented_nii = make_nii(segmented, [1 1 1],[8.3895 21.5912  25.1251] , 4); 
% save_nii(segmented_nii, [result_path_brain,'\',type,'_',name,'_segmented_tumor.nii.gz']);
%save_nii(segmented_nii, ['C:\Users\havm2701\Dropbox\Tumor_fortin\','segmented_tumor_FUSION_3.nii.gz']);


%save statistics
% save([result_path_brain,'\',type,'_',name,'_stat_results.mat'],'precision','recall');

%save meta files(.mha)
%mha_result_name = strrep(flair_name , 'XX.O.MR_Flair',['Seg_',type,'_',name]);
%mha_result_path = [result_path_brain,'\',mha_result_name];
 mha_result_path = ['HG_0001_withspatialinfo.mha'];
segmented=uint16(segmented);
 writemetaimagefile(mha_result_path, segmented, info.PixelDimensions,info.Offset);
% MASK = uint16(MASK);
% mha_MASK_path = [result_path_brain,'\','MASK_',type,'_',name,'.mha'];
% writemetaimagefile(mha_MASK_path, MASK, infot1c.PixelDimensions,infot1c.Offset);
