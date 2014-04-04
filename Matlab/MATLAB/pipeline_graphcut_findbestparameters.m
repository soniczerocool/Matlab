clear all
close all
mex -lut C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\interactiveGraphcut_withdataterms.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\GCoptimization.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\graph.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\LinkedBlockList.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\maxflow.cpp

stats = [];


global T1C T2 FLAIR
%% define paths and parameters
Nb_Selected = 500;

Nb_healthy = 500;
downsamplerate = 20;
 results_path_root = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\CHALLENGE2013_KNN_graphcut\';
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
 MASK = double(MASK);


 [ymin,ymax,xmin,xmax,zmin,zmax] = make_regtangle(T1C);
T1C = T1C(xmin:xmax , ymin:ymax , zmin:zmax);
T2 = T2(xmin:xmax , ymin:ymax , zmin:zmax);
FLAIR = FLAIR(xmin:xmax , ymin:ymax , zmin:zmax);

MASK = MASK(xmin:xmax , ymin:ymax , zmin:zmax);
truth_original = truth;
truth = truth(xmin:xmax , ymin:ymax , zmin:zmax);
%% make space and selected points space
[space,mask_idx,truth_idx] = make_space(T1C, FLAIR, T2,truth, MASK);% changed the way the space is created

[selected_space, mask_idx] = make_selected_space(space, mask_idx, Nb_healthy,Nb_Selected );


save_path = [results_path_root, 'statistics','_',name,'.txt'];


f = fopen(save_path,'w');

for k = 3:11
   
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
% space_test = space(1:2,:);
% space_test(1,:) = space(1,:)*max(T1C(:));
% space_test(2,:) = space(2,:)*max(T2(:));
    for alpha = 0.02:0.005:0.1
energy = alpha * energy;
 lf = interactiveGraphcut_withdataterms(space,MASK,energy,0.005);
 [h,w,d] = size(T1C);
 lf = reshape(lf,h, w, d);
 
  segmented_volume = lf ;
 segmented_brain = zeros(height,width,depth);
segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
segmented_brain(segmented_brain == 3)=4;
segmented_brain(segmented_brain == 1)=3;
truth = truth_original;
[dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented_brain, truth);
 %segmented_nii = make_nii(segmented); 
 %save_nii(segmented_nii, [result_path_brain,'\',type,'_',name,'_segmented_tumor.nii.gz']);
% 
% %save statistics
% save([result_path_brain,'\',type,'_',name,'_stat_results.mat'],'precision','recall');
% 
% %save meta files(.mha)
%mha_result_name = strrep(flair_name , 'XX.O.MR_Flair',['Seg_',type,'_',name]);
%mha_result_path = 'boykov_segment_4l_3_0.076abel.mha';
%segmented_brain=uint16(segmented_brain);
 %writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset);
 fprintf(f,['statistics for k_',num2str(k),'_alpha_',num2str(alpha)])
 fprintf(f,'\n');
 avg = write_txt_file_stats(precision,recall,dice,jaccard, f)
 stats = [stats , avg];
  fprintf(f,'\n');
 display(['k_',k,'_alpha',alpha])
    end
end
fclose(f)