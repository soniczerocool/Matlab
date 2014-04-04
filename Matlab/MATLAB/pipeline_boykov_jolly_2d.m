clear all
close all
mex -lut C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\boykov_jolly_2d_binary.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\GCoptimization.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\graph.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\LinkedBlockList.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\maxflow.cpp

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
MASK = make_mask(T1C, T2, FLAIR);

MASK = MASK(:,:,90);
% or load the mask from a path
 %load('MASK.mat');
MASK = double(MASK);
%% make space and selected points space

T1C = double(T1C(:,:,90));
T2 = double(T2(:,:,90));
FLAIR = double(FLAIR(:,:,90));
T1C = T1C/max(max(max(T1C)));
T2 = T2/max(max(max(T2)));
FLAIR = FLAIR/max(max(max(FLAIR)));


% 
% truth(truth<0.5)=0;                 %healthy 
% truth(0.5<truth & truth<1.5)=2;     %tomur
% truth(1.5<truth & truth<2.5)=1;     %edema

%% creating the 6 D space (T1,T2,FLAIR,x,y,z)

space     = zeros(5,length(T1C(:)));
mask_idx  = zeros(1,size(space,2));
truth_idx = zeros(1,size(space,2));
[height,width] = size(T1C);
i=1; 
for ROW = 1:height
   for COL=1:width
      
            
                space(:,i)=[T1C(ROW,COL);T2(ROW,COL);FLAIR(ROW,COL);ROW/height;COL/width];
                mask_idx(i) = MASK(ROW,COL);
                truth_idx(i) = truth(ROW,COL);
                i = i+1;
   end
end






%%
% background = find(sum(space(1:3,:))< 0.0001);
%  space(:,background) = [];
%   mask_idx(background) = [];
% 
% [selected_space, mask_idx] = make_selected_space(space, mask_idx, Nb_healthy,Nb_Selected );


 %% apply graph cut
 mask = MASK;
 image = T1C;
 lf = boykov_jolly_2d_binary(image,mask,60);
 [h,w] = size(T1C);
 lf = reshape(lf,h, w);

figure,imshow(lf,[])
figure, imshow(image,[])
figure, imshow(mask,[])
 %  
%   segmented_volume = lf ;
%  segmented_brain = zeros(height,width,depth);
% segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
% 
% 
%  %segmented_nii = make_nii(segmented); 
%  %save_nii(segmented_nii, [result_path_brain,'\',type,'_',name,'_segmented_tumor.nii.gz']);
% % 
% % %save statistics
% % save([result_path_brain,'\',type,'_',name,'_stat_results.mat'],'precision','recall');
% % 
% % %save meta files(.mha)
% %mha_result_name = strrep(flair_name , 'XX.O.MR_Flair',['Seg_',type,'_',name]);
% mha_result_path = 'boykov_segment.mha';
% segmented_brain=uint16(segmented_brain);
%  writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset);
%  