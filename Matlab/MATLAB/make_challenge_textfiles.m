clear all
close all;
addpath(genpath('C:\Users\havm2701\Dropbox\PhD\Matlab\MATLAB\MHA'))
run('C:\Users\havm2701\Dropbox\PhD\Matlab\MATLAB\vlfeat-0.9.17\toolbox\vl_setup');

global T1C T2 FLAIR 



k=0;
 %brain_dir ='N:\User\Research\Brain_data\BRATS-2\Image_Data\';
brain_dir = ['N:\User\Research\Brain_data\BRATS2013_CHALLENGE\Challenge\'];
  
%% define paths and parameters
Nb_Selected = 2000;

Nb_healthy = 2000;
downsamplerate = 20;
Nb_Brains = 30;
brain_type_list = dir(brain_dir);
brain_type_list = brain_type_list(3:end);
 for subfoldercounter=1:length(brain_type_list)
     type = brain_type_list(subfoldercounter).name;
     brain_type_path = [brain_dir,type];
     brain_list = dir(brain_type_path);
     brain_list =brain_list(3:end);
%      if type == 'HG'
%           continue
%       end
     for brain_counter = 1:length(brain_list)
         name = brain_list(brain_counter).name;
         brain_name_path = [brain_type_path,'\',name];
         

%% load modalities
        [T1C, T2, FLAIR, truth, info,flair_name] = load_modalities_mac(brain_name_path);
        [height width depth] = size(T1C);
        T1C = double(T1C);
        truth = zeros(size(T2,1),size(T2,2),size(T2,3));

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
 %mask_name_path = ['N:\User\Research\Brain_data\BRATS-2_MASKS\Image_Data\',type,'\',name];
 mask_name_path = ['N:\User\Research\Brain_data\BRATS2013_CHALLENGE_MASKS\Image_Data\',type,'\',name];
 MASK = load_mask(mask_name_path);
% you can Edit the mask
     %MASK = edit_mask(T1C, T2, FLAIR, MASK);
 
%[ymin,ymax,xmin,xmax,zmin,zmax] = make_regtangle(T1C);  
% save ROI (regtangle) coordinate points
 

%% load reg coordinates from the saved location
% roi_save_dir = ['N:\User\Research\Brain_data\BRATS-2_roi_coordinates\Image_Data\',type,'\',name];
%  roi_save_path = [roi_save_dir,'\ROI_',type,'_',name,'.mat'] ; 
% load(roi_save_path)
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
% 
% f = fopen(save_path1,'a');
%         
%         fprintf(f,'\n');
%         fprintf(f,'\n');
%         fprintf(f,'\n');
%         fprintf(f,'********************************************\n')
%         fprintf(f,['BRAIN_',type,'_',name]);
%        
%         fprintf(f,'\n');
%         fprintf(f,'*************')
%         fprintf(f,'\n');
%         f = fclose(f)




%%
% T1C = T1C(xmin:xmax , ymin:ymax , zmin:zmax);
% T2 = T2(xmin:xmax , ymin:ymax , zmin:zmax);
% FLAIR = FLAIR(xmin:xmax , ymin:ymax , zmin:zmax);
% 
% MASK = MASK(xmin:xmax , ymin:ymax , zmin:zmax);
% truth = truth(xmin:xmax , ymin:ymax , zmin:zmax);
 
% make space and selected space
 [space,mask_idx,truth_idx] = make_space(T1C, FLAIR, T2,truth, MASK);

[selected_space, mask_idx] = make_selected_space(space, mask_idx, Nb_healthy,Nb_Selected );
%%  make space and selected points space( MADE CHANGES FOR MAKING TEXT FILES >> FIX LATER)
  
[space,mask_idx,truth_idx]  = make_space(T1C, FLAIR, T2,truth, MASK);
         background = find(sum(space(1:3,:))< 0.001); 
         space(:,background) = [];
         truth_idx(background) = [];
          mask_idx(background) = [];
[selected_space, mask_idx] = make_selected_space(space, mask_idx, Nb_healthy,Nb_Selected );


 make_txt_file(space,selected_space,truth_idx,mask_idx,type, name);

        
%%  *****************COMMENTED FOR THE PURPOSE OF MAKE_TXT_FILE ********************************       
% % 
% % KNN search
% 
% kkk=k;
% [IDX,D] = knnsearch(selected_space',space','K', kkk ,'NSMethod','kdtree' ,'Distance','euclidean');
% 
% %[ind,knndist] = kdtreeKNN(selected_space,mask_idx,space,k);
% %IDX = ind';
% %D = knndist';
% 
%  [h,w,d] = size(T1C);
% % assign labels
% segmented_volume = assign_labels(IDX,D,mask_idx, space, h, w, d);
% 
%  segmented = zeros(height,width,depth);
% segmented(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
%    % compute statistics 
%  truth = truth_backup;
%  [dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented,truth);
%   stats = [precision;recall;dice;jaccard];
%   avg = sum(stats(:))/16;        
% 
% %%%% save mha file
% mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
% mha_result_path = [knn_result_save_dir,mha_result_name];
% segmented = uint16(segmented);
% writemetaimagefile(mha_result_path, segmented, info.PixelDimensions,info.Offset);
% %%%% save the matfile for statistics
% evaluation_matrix_savedir = [knn_result_save_dir,'evaluation_matrix\'];
% evaluation_matrix_savepath = [evaluation_matrix_savedir,type,'_',name,'.mat'];
% evaluation_matrix = [precision;recall;dice;jaccard];
% if exist(evaluation_matrix_savedir,'dir')==0
%     mkdir(evaluation_matrix_savedir)
% end
% save(evaluation_matrix_savepath,'evaluation_matrix','precision','recall','dice','jaccard')   
% %%%%%%%%%%% save results         
%         f = fopen(save_path1,'a');
%   
%         fprintf(f,['pure KNN k_',num2str(k)])
%         fprintf(f,'\n');
%         fprintf(f,'average statistics: ');
%         fprintf(f,num2str(avg));
%         fprintf(f,'\n');
%         fprintf(f,'*************')
%         fprintf(f,'\n');
% 
%         fclose(f);   
%         path_textfiles =textfiles_save_dir;
%  results_path_root_text = [path_textfiles,'\BRAIN_statistics\'];
%  if exist(results_path_root_text,'dir')==0
%     mkdir(results_path_root_text)
%  end
%  save_path2 = [results_path_root_text,type,'_',name,'_compare_methods.txt'];
%  g = fopen(save_path2,'a');
%             fprintf(f,'\n');
%         fprintf(f,'**************************************************')
%         fprintf(f,'\n');
%     fprintf(g,'pure KNN');
%     fprintf(g,['statistics for k_',num2str(kkk)]);
%      fprintf(g,'\n');
%      avg = write_txt_file_stats(precision,recall,dice,jaccard, g);
%     % stats = [stats , avg];
%       fprintf(g,'\n'); 
%       fclose(g);
% 
% %% apply median filter
%  segmented_volume = medfilt3(segmented_volume, [5,5,5]);
%  segmented = zeros(height,width,depth);
% segmented(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
%    % compute statistics 
%  truth = truth_backup;
%  [dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented,truth);
%   stats = [precision;recall;dice;jaccard];
%   avg = sum(stats(:))/16;        
% 
% %%%% save mha file
% mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
% mha_result_path = [knn_medianfilter_result_save_dir,mha_result_name];
% segmented = uint16(segmented);
% writemetaimagefile(mha_result_path, segmented, info.PixelDimensions,info.Offset);
% %%%% save the matfile for statistics
% evaluation_matrix_savedir = [knn_medianfilter_result_save_dir,'evaluation_matrix\'];
% evaluation_matrix_savepath = [evaluation_matrix_savedir,type,'_',name,'.mat'];
% evaluation_matrix = [precision;recall;dice;jaccard];
% if exist(evaluation_matrix_savedir,'dir')==0
%     mkdir(evaluation_matrix_savedir)
% end
% save(evaluation_matrix_savepath,'evaluation_matrix','precision','recall','dice','jaccard')   
% %%%%%%%%%%% save results         
%         f = fopen(save_path1,'a');
%   
%         fprintf(f,['KNN with median filter k_',num2str(k)])
%         fprintf(f,'\n');
%         fprintf(f,'average statistics: ');
%         fprintf(f,num2str(avg));
%         fprintf(f,'\n');
%         fprintf(f,'*************')
%         fprintf(f,'\n');
% 
%         fclose(f);   
%         path_textfiles =textfiles_save_dir;
%  results_path_root_text = [path_textfiles,'\BRAIN_statistics\'];
%  if exist(results_path_root_text,'dir')==0
%     mkdir(results_path_root_text)
%  end
%  save_path2 = [results_path_root_text,type,'_',name,'_compare_methods.txt'];
%  g = fopen(save_path2,'a');
%             fprintf(f,'\n');
%         fprintf(f,'***********************************************************************')
%         fprintf(f,'\n');
%     fprintf(g,'KNN with Median_filtering');
%     fprintf(g,['statistics for k_',num2str(kkk)]);
%      fprintf(g,'\n');
%      avg = write_txt_file_stats(precision,recall,dice,jaccard, g);
%     % stats = [stats , avg];
%       fprintf(g,'\n'); 
%       fclose(g);
%       
%  %% make data-term
% 
% matrix = make_dataterm_matrix(IDX,mask_idx,space);
% matrix = (matrix/k);
% %energy = -log(matrix);
% energy = 1 - matrix;
%  
%  
%  
% %% MRF method
% 
% 
% % T1C = T1C(xmin:xmax , ymin:ymax , zmin:zmax);
% % T2 = T2(xmin:xmax , ymin:ymax , zmin:zmax);
% % FLAIR = FLAIR(xmin:xmax , ymin:ymax , zmin:zmax);
% % 
% % MASK = MASK(xmin:xmax , ymin:ymax , zmin:zmax);
% % truth = truth(xmin:xmax , ymin:ymax , zmin:zmax);
% %  % make space and selected space
% %  [space,mask_idx,truth_idx] = make_space(T1C, FLAIR, T2,truth, MASK);
% % 
% % [selected_space, mask_idx] = make_selected_space(space, mask_idx, Nb_healthy,Nb_Selected );
% % % KNN search
% %[IDX,D] = knnsearch(selected_space',space','K', k ,'NSMethod','kdtree' ,'Distance','euclidean');
% 
% % make dataterm
% % matrix = make_dataterm_matrix(IDX,mask_idx,space);
% % matrix = (matrix/k);
% % %energy = -log(matrix);
% % energy = 1 - matrix;
% % apply graphcut MRF
% beta2 = 0.2;
% segmented = markov_network_multilable_general3Dgraph(double(T1C),energy,beta2);
% [h,w,d] = size(T1C);
% segmented_volume = label_to_matrix(segmented,space,h,w,d);
% % compute statistics 
% segmented_brain = zeros(height,width,depth);
% segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
% truth = truth_backup;
%  segmented_brain(segmented_brain == 3)=4;
%  segmented_brain(segmented_brain == 1)=3;
% [dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented_brain,truth);
%  stats = [precision;recall;dice;jaccard];
%         avg = sum(stats(:))/16;
% %%%%%%%%%%%%%% save results        
% mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
% mha_result_path = [mrf_result_save_dir,mha_result_name];
% segmented_brain = uint16(segmented_brain);
% writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset);        
% %%%% save the matfile for statistics
% evaluation_matrix_savedir = [mrf_result_save_dir,'evaluation_matrix\'];
% evaluation_matrix_savepath = [evaluation_matrix_savedir,type,'_',name,'.mat'];
% evaluation_matrix = [precision;recall;dice;jaccard];
% if exist(evaluation_matrix_savedir,'dir')==0
%     mkdir(evaluation_matrix_savedir)
% end
% save(evaluation_matrix_savepath,'evaluation_matrix','precision','recall','dice','jaccard')      
%          
%         f = fopen(save_path1,'a');
%         fprintf(f,['MRF without denoising denoising:'])
%    
%         fprintf(f,['( k_',num2str(k),'_beta2_',num2str(beta2),')'])
%         fprintf(f,'\n');
%         fprintf(f,'average_statistics: ');
%         fprintf(f,num2str(avg));
%         fprintf(f,'\n');
%         fprintf(f,'*************')
%         fprintf(f,'\n');
%         fclose(f);   
%      
%     g = fopen(save_path2,'a');
%     fprintf(g,'MRF_without denoising');
%     fprintf(g,['statistics for k_',num2str(k),'_beta2_',num2str(beta2)]);
%     fprintf(g,'\n');
%     avg = write_txt_file_stats(precision,recall,dice,jaccard, g);
%     fprintf(g,'\n');     
%     fclose(g);
% 
% % apply median filter
% segmented_volume = medfilt3(segmented_volume, [3,3,3]);
% % compute statistics 
% segmented_brain = zeros(height,width,depth);
% segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
% truth = truth_backup;
%  segmented_brain(segmented_brain == 3)=4;
%  segmented_brain(segmented_brain == 1)=3;
% [dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented_brain,truth);
%  stats = [precision;recall;dice;jaccard];
%         avg = sum(stats(:))/16;
% %%%%%%%%%% save results
% mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
% mha_result_path = [mrf_denoised_result_save_dir,mha_result_name];
% segmented_brain = uint16(segmented_brain);
% writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset);
% %%%% save the matfile for statistics
% evaluation_matrix_savedir = [mrf_denoised_result_save_dir,'evaluation_matrix\'];
% evaluation_matrix_savepath = [evaluation_matrix_savedir,type,'_',name,'.mat'];
% evaluation_matrix = [precision;recall;dice;jaccard];
% if exist(evaluation_matrix_savedir,'dir')==0
%     mkdir(evaluation_matrix_savedir)
% end
% save(evaluation_matrix_savepath,'evaluation_matrix','precision','recall','dice','jaccard')
% %%%%save text files
% f = fopen(save_path1,'a');
%         fprintf(f,['MRF after denoising[3,3,3]:'])
%         fprintf(f,'\n');
%         fprintf(f,[' k_',num2str(k),'_beta2_',num2str(beta2)])
%         fprintf(f,'\n');
%         fprintf(f,'average: ');
%         fprintf(f,num2str(avg));
%          fprintf(f,'\n');
%         fprintf(f,'*************')
%         fprintf(f,'\n');
% 
% fclose(f);   
%      
%     g = fopen(save_path2,'a');
%     fprintf(g,'MRF_after denoising');
%     fprintf(g,['statistics for k_',num2str(k),'_beta2_',num2str(beta2)]);
%     fprintf(g,'\n');
%     avg = write_txt_file_stats(precision,recall,dice,jaccard, g);
%     fprintf(g,'\n');     
%     fclose(g);
% %% APPLY graphcut BOYKOV JOLLY with dataterms
% alpha = 0.05;
%  energy = alpha * energy;
%  beta = .0051;
% lf = interactiveGraphcut_withdataterms(space,double(MASK),energy,beta);
%   [h,w,d] = size(T1C);
%   lf = reshape(lf,h, w, d);
%    segmented_volume = lf ;
% % compute statistics 
% segmented_brain = zeros(height,width,depth);
% segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
% 
%  segmented_brain(segmented_brain == 3)=4;
%  segmented_brain(segmented_brain == 1)=3;
% [dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented_brain,truth);
%  stats = [precision;recall;dice;jaccard];
%         avg = sum(stats(:))/16;
% %%%%%%%%%%%% save results
% mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
% mha_result_path = [boykov_result_save_dir,mha_result_name];
% segmented_brain = uint16(segmented_brain);
% writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset); 
% %%%% save the matfile for statistics
% evaluation_matrix_savedir = [boykov_result_save_dir,'evaluation_matrix\'];
% evaluation_matrix_savepath = [evaluation_matrix_savedir,type,'_',name,'.mat'];
% evaluation_matrix = [precision;recall;dice;jaccard];
% if exist(evaluation_matrix_savedir,'dir')==0
%     mkdir(evaluation_matrix_savedir)
% end
% save(evaluation_matrix_savepath,'evaluation_matrix','precision','recall','dice','jaccard')   
% 
% %%%%save textfiles
%         f = fopen(save_path1,'a');
%         fprintf(f,['boykov without denoising denoising:'])
%         fprintf(f,['(k_',num2str(k),'_alpha_',num2str(alpha),'_beta_',num2str(beta),')'])
%         fprintf(f,'average statistics: ');
%         fprintf(f,num2str(avg));
%         fprintf(f,'\n');
%         fprintf(f,'*************')
%         fprintf(f,'\n');
%        fclose(f);   
%      
%     g = fopen(save_path2,'a');
%     fprintf(g,'boykov_without denoising');
%     fprintf(g,['statistics for k_',num2str(k),'_alpha_',num2str(alpha),'_beta_',num2str(beta)]);
%     fprintf(g,'\n');
%     avg = write_txt_file_stats(precision,recall,dice,jaccard, g);
%     fprintf(g,'\n');     
%     fclose(g);
% 
% % apply median filter
% segmented_volume = medfilt3(segmented_volume, [3,3,3]);
% % compute statistics 
% segmented_brain = zeros(height,width,depth);
% segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
% truth = truth_backup;
%  segmented_brain(segmented_brain == 3)=4;
%  segmented_brain(segmented_brain == 1)=3;
% [dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented_brain,truth);
%  stats = [precision;recall;dice;jaccard];
%         avg = sum(stats(:))/16;
% %%%%%%%% save results 
% mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
% mha_result_path = [boykov_denoised_result_save_dir,mha_result_name];
% segmented_brain = uint16(segmented_brain);
% writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset); 
% %%%% save the matfile for statistics
% evaluation_matrix_savedir = [boykov_denoised_result_save_dir,'evaluation_matrix\'];
% evaluation_matrix_savepath = [evaluation_matrix_savedir,type,'_',name,'.mat'];
% evaluation_matrix = [precision;recall;dice;jaccard];
% if exist(evaluation_matrix_savedir,'dir')==0
%     mkdir(evaluation_matrix_savedir)
% end
% save(evaluation_matrix_savepath,'evaluation_matrix','precision','recall','dice','jaccard')   
% 
% 
%         f = fopen(save_path1,'a');
%         fprintf(f,['boykov after denoising[3,3,3]:']) 
%         fprintf(f,['( k_',num2str(k),'_alpha_',num2str(alpha),'_beta_',num2str(beta),')'])
%         fprintf(f,'\n');
%         fprintf(f,'average_Statistic: ');
%         fprintf(f,num2str(avg));
%         fprintf(f,'\n');
%         fprintf(f,'*************')
%         fprintf(f,'\n');
%         fclose(f);   
%      
%     g = fopen(save_path2,'a');
%     fprintf(g,'boykov_after denoising');
%     fprintf(g,['statistics for k_',num2str(k),'_alpha_',num2str(alpha),'_beta_',num2str(beta)]);
%     fprintf(g,'\n');
%     avg = write_txt_file_stats(precision,recall,dice,jaccard, g);
%     fprintf(g,'\n');     
%     fclose(g);
%     
%     %% APPLY graphcut BOYKOV JOLLY without dataterms
% % 
%  energy2 = 0 * energy;
%  beta = .0051;
% lf = interactiveGraphcut_withdataterms(space,double(MASK),energy2,beta);
%   [h,w,d] = size(T1C);
%   lf = reshape(lf,h, w, d);
%    segmented_volume = lf ;
% % compute statistics 
% segmented_brain = zeros(height,width,depth);
% segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
% truth = truth_backup;
%  segmented_brain(segmented_brain == 3)=4;
%  segmented_brain(segmented_brain == 1)=3;
% [dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented_brain,truth);
%  stats = [precision;recall;dice;jaccard];
%         avg = sum(stats(:))/16;
% %%%%%%%%%%% save results
% mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
% mha_result_path = [bykov_ndt_result_save_dir,mha_result_name];
% segmented_brain = uint16(segmented_brain);
% writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset); 
% %%%% save the matfile for statistics
% evaluation_matrix_savedir = [bykov_ndt_result_save_dir,'evaluation_matrix\'];
% evaluation_matrix_savepath = [evaluation_matrix_savedir,type,'_',name,'.mat'];
% evaluation_matrix = [precision;recall;dice;jaccard];
% if exist(evaluation_matrix_savedir,'dir')==0
%     mkdir(evaluation_matrix_savedir)
% end
% save(evaluation_matrix_savepath,'evaluation_matrix','precision','recall','dice','jaccard')   
% 
%         f = fopen(save_path1,'a');
%         fprintf(f,['boykov(no dataterm) without denoising denoising:'])
%         fprintf(f,['(k_',num2str(k),'_alpha_',num2str(alpha),'_beta_',num2str(beta),')'])
%         fprintf(f,'average statistics: ');
%         fprintf(f,num2str(avg));
%         fprintf(f,'\n');
%         fprintf(f,'*************')
%         fprintf(f,'\n');
%        fclose(f);   
%      
%     g = fopen(save_path2,'a');
%     fprintf(g,'boykov(no dataterm)_without denoising');
%     fprintf(g,['statistics for k_',num2str(k),'_alpha_',num2str(alpha),'_beta_',num2str(beta)]);
%     fprintf(g,'\n');
%     avg = write_txt_file_stats(precision,recall,dice,jaccard, g);
%     fprintf(g,'\n');     
%     fclose(g);
% 
% % apply median filter
% segmented_volume = medfilt3(segmented_volume, [3,3,3]);
% % compute statistics 
% segmented_brain = zeros(height,width,depth);
% segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
% truth = truth_backup;
%  segmented_brain(segmented_brain == 3)=4;
%  segmented_brain(segmented_brain == 1)=3;
% [dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented_brain,truth);
%  stats = [precision;recall;dice;jaccard];
%         avg = sum(stats(:))/16;
% %%%%%%%%%%%%% save results
% mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
% mha_result_path = [bykov_ndt_denoised_result_save_dir,mha_result_name];
% segmented_brain = uint16(segmented_brain);
% writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset); 
% %%%% save the matfile for statistics
% evaluation_matrix_savedir = [bykov_ndt_denoised_result_save_dir,'evaluation_matrix\'];
% evaluation_matrix_savepath = [evaluation_matrix_savedir,type,'_',name,'.mat'];
% evaluation_matrix = [precision;recall;dice;jaccard];
% if exist(evaluation_matrix_savedir,'dir')==0
%     mkdir(evaluation_matrix_savedir)
% end
% save(evaluation_matrix_savepath,'evaluation_matrix','precision','recall','dice','jaccard')   
% 
% 
%         f = fopen(save_path1,'a');
%         fprintf(f,['boykov(no dataterm) after denoising[3,3,3]:']) 
%         fprintf(f,['( k_',num2str(k),'_alpha_',num2str(alpha),'_beta_',num2str(beta),')'])
%         fprintf(f,'\n');
%         fprintf(f,'average_Statistic: ');
%         fprintf(f,num2str(avg));
%         fprintf(f,'\n');
%         fprintf(f,'*************')
%         fprintf(f,'\n');
%         fclose(f);   
%      
%     g = fopen(save_path2,'a');
%     fprintf(g,'boykov(no dataterm)_after denoising');
%     fprintf(g,['statistics for k_',num2str(k),'_alpha_',num2str(alpha),'_beta_',num2str(beta)]);
%     fprintf(g,'\n');
%     avg = write_txt_file_stats(precision,recall,dice,jaccard, g);
%     fprintf(g,'\n');     
%     fclose(g);
  
     end
 end