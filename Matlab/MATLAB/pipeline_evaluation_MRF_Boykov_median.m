clear all
close all;
mex -lut C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\markov_network_multilable_general3Dgraph.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\GCoptimization.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\graph.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\LinkedBlockList.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\maxflow.cpp
mex -lut C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\interactiveGraphcut_withdataterms.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\GCoptimization.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\graph.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\LinkedBlockList.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\maxflow.cpp


Nb_Selected = 500;
Nb_healthy = 500;
%result_path_brain = 'J:\User\Research\Results\BRATS-2_KNN_ENDAUGUST_FORUPLOAD\';
brain_dir ='N:\User\Research\Brain_data\BRATS-2\Image_Data\';
%brain_dir = '/Volumes/Expansion Drive/User/Research/Brain_data/BRATS-2/Image_Data/';
%result_path_brain = 'N:\User\Research\Results\textfiles_oct\';
%mask_dir = '/Users/uoft/Dropbox/PhD/BratsAnalysis/CHALLENGE2013_KNN_mha/';
mask_dir = 'N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\';
result_path = 'N:\User\Research\Results\BRATS-2_roi_coordinates\';
path_textfiles = 'N:\User\Research\Results\BRATS-2_endNov_textresult\';
%result_path = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\KNN_stat_result_for_k\';
% result_path = 'J:\User\Research\Results\BRATS-2_KNN_ENDAUGUST_lowresT2\';
if exist(result_path,'dir')==0
    mkdir(result_path)
end
if exist(path_textfiles,'dir')==0
    mkdir(path_textfiles)
end

save_path1 = [path_textfiles,'summery_evaluation_methods.txt'];
% if exist(result_path_brain,'dir')==0
%     mkdir(result_path_brain)
% end

Nb_Brains = 30;
mask_type_list = dir(mask_dir);
mask_type_list = mask_type_list(3:end);
  
 
 precision_total = zeros(1,30);
 recall_total = zeros(1,30);
 dice_total = zeros(1,30);
 
 for kkk= 3:3
 INDEX = 1;
 if kkk==0
     kkk=1;
 end
 
 precesion_recall = zeros(2,Nb_Brains);
 for subfoldercounter=1:length(mask_type_list)
     type = mask_type_list(subfoldercounter).name;
     brain_type_path = [brain_dir,type];
     brain_list = dir(brain_type_path);
     brain_list =brain_list(3:end);
     mask_type_path = [mask_dir,type];
     mask_list = dir(mask_type_path);
     mask_list = mask_list(4:end);
     for brain_counter = 12:length(mask_list)
         name = brain_list(brain_counter).name;
         name2 = mask_list(brain_counter).name;
         brain_name_path = [brain_type_path,'\',name];
         mask_name_path = [mask_type_path,'\',name];
         TF = strcmp(name,name2);
         if TF == 0
             error('Brain name does not match the mask name!')
         end
         %% start processing for each brain with id [type,'_',name]
         t_start = tic;
        disp(['Processing brain ',type,'_',name, ' with k = ', num2str(kkk)])
        f = fopen(save_path1,'a');
        fprintf(f,['**************************************************']);
        fprintf(f,'\n')
        fprintf(f,['Brain _',type,'_',name]);
        fprintf(f,'\n');
        fprintf(f,'*************')
        fprintf(f,'\n');
        fclose(f);
    
%% load modalities
        [T1C, T2, FLAIR, truth, info,flair_name] = load_modalities_mac(brain_name_path);
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
        
% creating the mask
        %MASK = make_mask(T1C, T2, FLAIR);
        % or load the mask from a path 
         MASK = load_mask(mask_name_path);
         

% make space and selected points space
       [space,mask_idx,truth_idx]  = make_space(T1C, FLAIR, T2,truth, MASK);
        

        %make_txt_file(space,truth_idx,mask_idx,type, name);
         background = find(sum(space(1:3,:))< 0.001); 
         space(:,background) = [];
         mask_idx(background) = [];
        
         [selected_space, mask_idx] = make_selected_space(space, mask_idx, Nb_healthy,Nb_Selected );
% 
% KNN search

k=kkk;
[IDX,D] = knnsearch(selected_space',space','K', kkk ,'NSMethod','kdtree' ,'Distance','euclidean');

% assign labels
segmented = assign_labels(IDX,D,mask_idx, space, height, width, depth);
% apply median filter
        segmented = medfilt3(segmented, [5,5,5]);
   % compute statistics 
        [dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented,truth);
        stats = [precision;recall;dice;jaccard];
        avg = sum(stats(:))/16;
% save results         
        f = fopen(save_path1,'a');
        fprintf(f,['Median k_',num2str(k)])
        fprintf(f,'\n');
        fprintf(f,'average statistics: ');
        fprintf(f,num2str(avg));
        fprintf(f,'\n');
        fprintf(f,'*************')
        fprintf(f,'\n');

        fclose(f);   
 results_path_root_text = [path_textfiles,type,'\',name];
 if exist(results_path_root_text,'dir')==0
    mkdir(results_path_root_text)
 end
 save_path2 = [results_path_root_text,'\compare_methods.txt'];
 g = fopen(save_path2,'a');
    fprintf(g,'Median_filtering');
    fprintf(g,['statistics for k_',num2str(kkk)]);
     fprintf(g,'\n');
     avg = write_txt_file_stats(precision,recall,dice,jaccard, g);
    % stats = [stats , avg];
      fprintf(g,'\n'); 
      fclose(g);
%% MRF method

 [ymin,ymax,xmin,xmax,zmin,zmax] = make_regtangle(T1C);
T1C = T1C(xmin:xmax , ymin:ymax , zmin:zmax);
T2 = T2(xmin:xmax , ymin:ymax , zmin:zmax);
FLAIR = FLAIR(xmin:xmax , ymin:ymax , zmin:zmax);

MASK = MASK(xmin:xmax , ymin:ymax , zmin:zmax);
truth = truth(xmin:xmax , ymin:ymax , zmin:zmax);
 % make space and selected space
 [space,mask_idx,truth_idx] = make_space(T1C, FLAIR, T2,truth, MASK);

[selected_space, mask_idx] = make_selected_space(space, mask_idx, Nb_healthy,Nb_Selected );
% KNN search
[IDX,D] = knnsearch(selected_space',space','K', k ,'NSMethod','kdtree' ,'Distance','euclidean');

% make dataterm
matrix = make_dataterm_matrix(IDX,mask_idx,space);
matrix = (matrix/k);
%energy = -log(matrix);
energy = 1 - matrix;
% apply graphcut MRF
beta2 = 0.2;
segmented = markov_network_multilable_general3Dgraph(T1C,energy,beta2);
[h,w,d] = size(T1C);
segmented_volume = label_to_matrix(segmented,space,h,w,d);
% compute statistics 
segmented_brain = zeros(height,width,depth);
segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
truth = truth_backup;
 segmented_brain(segmented_brain == 3)=4;
 segmented_brain(segmented_brain == 1)=3;
[dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented_brain,truth);
 stats = [precision;recall;dice;jaccard];
        avg = sum(stats(:))/16;
% save results         
        f = fopen(save_path1,'a');
        fprintf(f,['MRF without denoising denoising:'])
   
        fprintf(f,['( k_',num2str(k),'_beta2_',num2str(beta2),')'])
        fprintf(f,'\n');
        fprintf(f,'average_statistics: ');
        fprintf(f,num2str(avg));
        fprintf(f,'\n');
        fprintf(f,'*************')
        fprintf(f,'\n');
        fclose(f);   
     
    g = fopen(save_path2,'a');
    fprintf(g,'MRF_without denoising');
    fprintf(g,['statistics for k_',num2str(k),'_beta2_',num2str(beta2)]);
    fprintf(g,'\n');
    avg = write_txt_file_stats(precision,recall,dice,jaccard, g);
    fprintf(g,'\n');     
    fclose(g);

% apply median filter
segmented_volume = medfilt3(segmented_volume, [3,3,3]);
% compute statistics 
segmented_brain = zeros(height,width,depth);
segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
truth = truth_backup;
 segmented_brain(segmented_brain == 3)=4;
 segmented_brain(segmented_brain == 1)=3;
[dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented_brain,truth);
 stats = [precision;recall;dice;jaccard];
        avg = sum(stats(:))/16;
% save results         
        f = fopen(save_path1,'a');
        fprintf(f,['MRF after denoising[3,3,3]:'])
        fprintf(f,'\n');
        fprintf(f,[' k_',num2str(k),'_beta2_',num2str(beta2)])
        fprintf(f,'\n');
        fprintf(f,'average: ');
        fprintf(f,num2str(avg));
         fprintf(f,'\n');
        fprintf(f,'*************')
        fprintf(f,'\n');

        fclose(f);   
     
    g = fopen(save_path2,'a');
    fprintf(g,'MRF_after denoising');
    fprintf(g,['statistics for k_',num2str(k),'_beta2_',num2str(beta2)]);
    fprintf(g,'\n');
    avg = write_txt_file_stats(precision,recall,dice,jaccard, g);
    fprintf(g,'\n');     
    fclose(g);
%% APPLY graphcut BOYKOV JOLLY with dataterms
alpha = 0.01;
 energy = alpha * energy;
 beta = .0051;
lf = interactiveGraphcut_withdataterms(space,double(MASK),energy,beta);
  [h,w,d] = size(T1C);
  lf = reshape(lf,h, w, d);
   segmented_volume = lf ;
% compute statistics 
segmented_brain = zeros(height,width,depth);
segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
truth = truth_backup;
 segmented_brain(segmented_brain == 3)=4;
 segmented_brain(segmented_brain == 1)=3;
[dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented_brain,truth);
 stats = [precision;recall;dice;jaccard];
        avg = sum(stats(:))/16;
% save results         
        f = fopen(save_path1,'a');
        fprintf(f,['boykov without denoising denoising:'])
        fprintf(f,['(k_',num2str(k),'_alpha_',num2str(alpha),'_beta_',num2str(beta),')'])
        fprintf(f,'average statistics: ');
        fprintf(f,num2str(avg));
        fprintf(f,'\n');
        fprintf(f,'*************')
        fprintf(f,'\n');
       fclose(f);   
     
    g = fopen(save_path2,'a');
    fprintf(g,'boykov_without denoising');
    fprintf(g,['statistics for k_',num2str(k),'_alpha_',num2str(alpha),'_beta_',num2str(beta)]);
    fprintf(g,'\n');
    avg = write_txt_file_stats(precision,recall,dice,jaccard, g);
    fprintf(g,'\n');     
    fclose(g);

% apply median filter
segmented_volume = medfilt3(segmented_volume, [3,3,3]);
% compute statistics 
segmented_brain = zeros(height,width,depth);
segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
truth = truth_backup;
 segmented_brain(segmented_brain == 3)=4;
 segmented_brain(segmented_brain == 1)=3;
[dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented_brain,truth);
 stats = [precision;recall;dice;jaccard];
        avg = sum(stats(:))/16;
% save results         
        f = fopen(save_path1,'a');
        fprintf(f,['boykov after denoising[3,3,3]:']) 
        fprintf(f,['( k_',num2str(k),'_alpha_',num2str(alpha),'_beta_',num2str(beta),')'])
        fprintf(f,'\n');
        fprintf(f,'average_Statistic: ');
        fprintf(f,num2str(avg));
        fprintf(f,'\n');
        fprintf(f,'*************')
        fprintf(f,'\n');
        fclose(f);   
     
    g = fopen(save_path2,'a');
    fprintf(g,'boykov_after denoising');
    fprintf(g,['statistics for k_',num2str(k),'_alpha_',num2str(alpha),'_beta_',num2str(beta)]);
    fprintf(g,'\n');
    avg = write_txt_file_stats(precision,recall,dice,jaccard, g);
    fprintf(g,'\n');     
    fclose(g);
    
    %% APPLY graphcut BOYKOV JOLLY without dataterms

 beta = .0051;
lf = interactiveGraphcut_withdataterms(space,double(MASK),beta);
  [h,w,d] = size(T1C);
  lf = reshape(lf,h, w, d);
   segmented_volume = lf ;
% compute statistics 
segmented_brain = zeros(height,width,depth);
segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
truth = truth_backup;
 segmented_brain(segmented_brain == 3)=4;
 segmented_brain(segmented_brain == 1)=3;
[dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented_brain,truth);
 stats = [precision;recall;dice;jaccard];
        avg = sum(stats(:))/16;
% save results         
        f = fopen(save_path1,'a');
        fprintf(f,['boykov(no dataterm) without denoising denoising:'])
        fprintf(f,['(k_',num2str(k),'_alpha_',num2str(alpha),'_beta_',num2str(beta),')'])
        fprintf(f,'average statistics: ');
        fprintf(f,num2str(avg));
        fprintf(f,'\n');
        fprintf(f,'*************')
        fprintf(f,'\n');
       fclose(f);   
     
    g = fopen(save_path2,'a');
    fprintf(g,'boykov(no data term)_without denoising');
    fprintf(g,['statistics for k_',num2str(k),'_alpha_',num2str(alpha),'_beta_',num2str(beta)]);
    fprintf(g,'\n');
    avg = write_txt_file_stats(precision,recall,dice,jaccard, g);
    fprintf(g,'\n');     
    fclose(g);

% apply median filter
segmented_volume = medfilt3(segmented_volume, [3,3,3]);
% compute statistics 
segmented_brain = zeros(height,width,depth);
segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
truth = truth_backup;
 segmented_brain(segmented_brain == 3)=4;
 segmented_brain(segmented_brain == 1)=3;
[dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented_brain,truth);
 stats = [precision;recall;dice;jaccard];
        avg = sum(stats(:))/16;
% save results         
        f = fopen(save_path1,'a');
        fprintf(f,['boykov(no dataterm) after denoising[3,3,3]:']) 
        fprintf(f,['( k_',num2str(k),'_alpha_',num2str(alpha),'_beta_',num2str(beta),')'])
        fprintf(f,'\n');
        fprintf(f,'average_Statistic: ');
        fprintf(f,num2str(avg));
        fprintf(f,'\n');
        fprintf(f,'*************')
        fprintf(f,'\n');
        fclose(f);   
     
    g = fopen(save_path2,'a');
    fprintf(g,'boykov(no dataterm)_after denoising');
    fprintf(g,['statistics for k_',num2str(k),'_alpha_',num2str(alpha),'_beta_',num2str(beta)]);
    fprintf(g,'\n');
    avg = write_txt_file_stats(precision,recall,dice,jaccard, g);
    fprintf(g,'\n');     
    fclose(g);
  %% save ROI (regtangle) coordinate points
  
 results_path_root = [result_path,type,'\',name];
 if exist(results_path_root,'dir')==0
    mkdir(results_path_root)
 end
   save([results_path_root,'\','coordinate_points_for_',type,'_',name,'.mat'],'ymin','ymax','xmin','xmax','zmin','zmax');



 end
 %% save statistical results
 %save([result_path,'stat_for_k_',num2str(kkk),'.mat'],'precesion_recall');
 
 end
 end