clear all
close all;
addpath(genpath('C:\Users\havm2701\Dropbox\PhD\Matlab\MATLAB\MHA'))
run('C:\Users\havm2701\Dropbox\PhD\Matlab\MATLAB\vlfeat-0.9.17\toolbox\vl_setup');
mex -lut C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\markov_network_multilable_general3Dgraph.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\GCoptimization.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\graph.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\LinkedBlockList.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\maxflow.cpp
mex -lut C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\interactiveGraphcut_withdataterms.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\GCoptimization.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\graph.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\LinkedBlockList.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\maxflow.cpp
global T1C T2 FLAIR 

k=3;

% output directories
main_result_dir = 'N:\User\Research\Results\BRATS-2_ICPR_BRATS_CHALLENGEDATA\';
textfiles_save_dir =[main_result_dir,'textresult\'];


%brain_dir = '/Volumes/Expansion Drive/User/Research/Brain_data/BRATS-2/Image_Data/';
%result_path_brain = 'N:\User\Research\Results\textfiles_oct\';
%mask_dir = '/Users/uoft/Dropbox/PhD/BratsAnalysis/CHALLENGE2013_KNN_mha/';


%result_path = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\KNN_stat_result_for_k\';
% result_path = 'J:\User\Research\Results\BRATS-2_KNN_ENDAUGUST_lowresT2\';


if exist(textfiles_save_dir,'dir')==0
    mkdir(textfiles_save_dir)
end

save_path1 = [textfiles_save_dir,'summery_evaluation_methods.txt'];
% if exist(result_path_brain,'dir')==0
%     mkdir(result_path_brain)
% end

%% define paths and parameters
Nb_Selected = 2000;

Nb_healthy = 2000;
downsamplerate = 20;


knn_result_save_dir =[main_result_dir,'knn\'];
if exist(knn_result_save_dir,'dir')==0
    mkdir(knn_result_save_dir)
end


knn_medianfilter_result_save_dir =[main_result_dir,'knn_medianfilter\'];
if exist(knn_medianfilter_result_save_dir,'dir')==0
    mkdir(knn_medianfilter_result_save_dir)
end

mrf_result_save_dir =  [main_result_dir,'mrf\nodenoising\'];
if exist(mrf_result_save_dir,'dir')==0
    mkdir(mrf_result_save_dir)
end

boykov_result_save_dir =  [main_result_dir,'boykov\nodenoising\'];
if exist(boykov_result_save_dir,'dir')==0
    mkdir(boykov_result_save_dir)
end

bykov_ndt_result_save_dir =  [main_result_dir,'boykov_nodataterm\nodenoising\'];
if exist(bykov_ndt_result_save_dir,'dir')==0
    mkdir(bykov_ndt_result_save_dir)
end

mrf_denoised_result_save_dir =  [main_result_dir,'mrf\Denoised\'];
if exist(mrf_denoised_result_save_dir,'dir')==0
    mkdir(mrf_denoised_result_save_dir)
end

boykov_denoised_result_save_dir =  [main_result_dir,'boykov\Denoised\'];
if exist(boykov_denoised_result_save_dir,'dir')==0
    mkdir(boykov_denoised_result_save_dir)
end

bykov_ndt_denoised_result_save_dir =  [main_result_dir,'boykov_nodataterm\Denoised\'];
if exist(bykov_ndt_denoised_result_save_dir,'dir')==0
    mkdir(bykov_ndt_denoised_result_save_dir)
end


%  brain_name_path = ['N:\User\Research\Brain_data\BRATS-2\Image_Data\',type,'\',name];
%  mask_name_path = ['N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\',type,'\',name];

 brain_dir = ['N:\User\Research\Brain_data\BRATS2013_CHALLENGE\Challenge\'];

 

Nb_Brains = 30;
brain_type_list = dir(brain_dir);
brain_type_list = brain_type_list(3:end);
 for subfoldercounter=1:length(brain_type_list)
     type = brain_type_list(subfoldercounter).name;
     brain_type_path = [brain_dir,type];
     brain_list = dir(brain_type_path);
     brain_list =brain_list(3:end);
    
     for brain_counter = 1:length(brain_list)
         name = brain_list(brain_counter).name;
         brain_name_path = [brain_type_path,'\',name];
         

%% load modalities
        [T1C, T2, FLAIR, truth, info,flair_name] = load_modalities_mac(brain_name_path);
        truth = zeros(size(T2,1),size(T2,2),size(T2,3));
        [height width depth] = size(T1C);
        T1C = double(T1C);
        T2 = double(T2);
        FLAIR = double(FLAIR);
        % take backup of the modalities
        truth_backup = truth;
        T1C_backup = T1C;
        T2_backup = T2;
        FLAIR_backup = FLAIR;
%% reduce the resolution of T2 (check for Maxime)
%        T2 = rescalingT2(T2);
        
%% creating the mask and ROI
   % MASK = make_mask(T1C, T2, FLAIR);
% or load the mask from a path
%mask_name_path = 'N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\HG\0026';
 mask_name_path = ['N:\User\Research\Brain_data\BRATS2013_CHALLENGE_MASKS\Image_Data\',type,'\',name];
 MASK = load_mask(mask_name_path);
% you can Edit the mask
     %MASK = edit_mask(T1C, T2, FLAIR, MASK);
 
%[ymin,ymax,xmin,xmax,zmin,zmax] = make_regtangle(T1C);  
% save ROI (regtangle) coordinate points
 

%% load reg coordinates from the saved location
roi_save_dir = ['N:\User\Research\Brain_data\BRATS2013_CHALLENGE_roi_coordinates\Image_Data\',type,'\',name];
 roi_save_path = [roi_save_dir,'\ROI_',type,'_',name,'.mat'] ; 
load(roi_save_path)
%% 
% %  if exist(roi_save_dir,'dir')==0
%     mkdir(roi_save_dir)
%  end
%    save(roi_save_path,'ymin','ymax','xmin','xmax','zmin','zmax');

 
% mask_save_path = [mask_save_dir,'\MASK_',type,'_',name,'.mha'] ;     
% MASK = uint16(MASK);
% writemetaimagefile(mask_save_path, MASK, info.PixelDimensions,info.Offset);     
% MASK = double(MASK);  
disp('***********************************')
disp(['Processing brain ',type,'_',name, ' with k = ', num2str(k)])
disp('***********************************')

%% write in the textfile

f = fopen(save_path1,'a');
        
        fprintf(f,'\n');
        fprintf(f,'\n');
        fprintf(f,'\n');
        fprintf(f,'********************************************\n')
        fprintf(f,['BRAIN_',type,'_',name]);
       
        fprintf(f,'\n');
        fprintf(f,'*************')
        fprintf(f,'\n');
        f = fclose(f)




%%
T1C = T1C(xmin:xmax , ymin:ymax , zmin:zmax);
T2 = T2(xmin:xmax , ymin:ymax , zmin:zmax);
FLAIR = FLAIR(xmin:xmax , ymin:ymax , zmin:zmax);

MASK = MASK(xmin:xmax , ymin:ymax , zmin:zmax);
truth = truth(xmin:xmax , ymin:ymax , zmin:zmax);
 
% make space and selected space
 [space,mask_idx,truth_idx] = make_space(T1C, FLAIR, T2,truth, MASK);

[selected_space, mask_idx] = make_selected_space(space, mask_idx, Nb_healthy,Nb_Selected );

  


        
%%  *****************COMMENTED FOR THE PURPOSE OF MAKE_TXT_FILE ********************************       
% 
% KNN search

kkk=k;
[IDX,D] = knnsearch(selected_space',space','K', kkk ,'NSMethod','kdtree' ,'Distance','euclidean');

%[ind,knndist] = kdtreeKNN(selected_space,mask_idx,space,k);
%IDX = ind';
%D = knndist';

 [h,w,d] = size(T1C);
% assign labels
segmented_volume = assign_labels(IDX,D,mask_idx, space, h, w, d);

 segmented = zeros(height,width,depth);
segmented(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
   % compute statistics 
 

%%%% save mha file
mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
mha_result_path = [knn_result_save_dir,mha_result_name];
segmented = uint16(segmented);
writemetaimagefile(mha_result_path, segmented, info.PixelDimensions,info.Offset);


%% apply median filter
 segmented_volume = medfilt3(segmented_volume, [5,5,5]);
 segmented = zeros(height,width,depth);
segmented(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
   
      

%%%% save mha file
mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
mha_result_path = [knn_medianfilter_result_save_dir,mha_result_name];
segmented = uint16(segmented);
writemetaimagefile(mha_result_path, segmented, info.PixelDimensions,info.Offset);

      



      
 %% make data-term

matrix = make_dataterm_matrix(IDX,mask_idx,space);
matrix = (matrix/k);
%energy = -log(matrix);
energy = 1 - matrix;
 
 
 
%% MRF method



beta2 = 0.2;
segmented = markov_network_multilable_general3Dgraph(double(T1C),energy,beta2);
[h,w,d] = size(T1C);
segmented_volume = label_to_matrix(segmented,space,h,w,d);
% compute statistics 
segmented_brain = zeros(height,width,depth);
segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
truth = truth_backup;
 segmented_brain(segmented_brain == 3)=4;
 segmented_brain(segmented_brain == 1)=3;

%%%%%%%%%%%%%% save results        
mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
mha_result_path = [mrf_result_save_dir,mha_result_name];
segmented_brain = uint16(segmented_brain);
writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset);        


% apply median filter
segmented_volume = medfilt3(segmented_volume, [3,3,3]);
% compute statistics 
segmented_brain = zeros(height,width,depth);
segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
truth = truth_backup;
 segmented_brain(segmented_brain == 3)=4;
 segmented_brain(segmented_brain == 1)=3;

%%%%%%%%%% save results
mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
mha_result_path = [mrf_denoised_result_save_dir,mha_result_name];
segmented_brain = uint16(segmented_brain);
writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset);

%% APPLY graphcut BOYKOV JOLLY with dataterms
alpha = 0.05;
 energy = alpha * energy;
 beta = .0051;
lf = interactiveGraphcut_withdataterms(space,double(MASK),energy,beta);
  [h,w,d] = size(T1C);
  lf = reshape(lf,h, w, d);
   segmented_volume = lf ;
% compute statistics 
segmented_brain = zeros(height,width,depth);
segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;

 segmented_brain(segmented_brain == 3)=4;
 segmented_brain(segmented_brain == 1)=3;

%%%%%%%%%%%% save results
mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
mha_result_path = [boykov_result_save_dir,mha_result_name];
segmented_brain = uint16(segmented_brain);
writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset); 



% apply median filter
segmented_volume = medfilt3(segmented_volume, [3,3,3]);
% compute statistics 
segmented_brain = zeros(height,width,depth);
segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
truth = truth_backup;
 segmented_brain(segmented_brain == 3)=4;
 segmented_brain(segmented_brain == 1)=3;

%%%%%%%% save results 
mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
mha_result_path = [boykov_denoised_result_save_dir,mha_result_name];
segmented_brain = uint16(segmented_brain);
writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset); 
  



    
    %% APPLY graphcut BOYKOV JOLLY without dataterms
% 
 energy2 = 0 * energy;
 beta = .0051;
lf = interactiveGraphcut_withdataterms(space,double(MASK),energy2,beta);
  [h,w,d] = size(T1C);
  lf = reshape(lf,h, w, d);
   segmented_volume = lf ;
% compute statistics 
segmented_brain = zeros(height,width,depth);
segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;

 segmented_brain(segmented_brain == 3)=4;
 segmented_brain(segmented_brain == 1)=3;

%%%%%%%%%%% save results
mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
mha_result_path = [bykov_ndt_result_save_dir,mha_result_name];
segmented_brain = uint16(segmented_brain);
writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset); 


 

% apply median filter
segmented_volume = medfilt3(segmented_volume, [3,3,3]);
% compute statistics 
segmented_brain = zeros(height,width,depth);
segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
truth = truth_backup;
 segmented_brain(segmented_brain == 3)=4;
 segmented_brain(segmented_brain == 1)=3;

 
%%%%%%%%%%%%% save results
mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
mha_result_path = [bykov_ndt_denoised_result_save_dir,mha_result_name];
segmented_brain = uint16(segmented_brain);
writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset); 
%%%% save the matfile for statistics

     end
 end