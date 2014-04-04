clear all
W=3;
Nb_Selected = 500;
Nb_healthy = 500;
%result_path_brain = 'J:\User\Research\Results\BRATS-2_KNN_ENDAUGUST_FORUPLOAD\';
result_path_brain = 'N:\User\Research\Results\BRATS-2_KNN_endDEC29test\';
if exist(result_path_brain,'dir')==0
    mkdir(result_path_brain)
end
brain_dir ='N:\User\Research\Brain_data\BRATS-2\Image_Data\';
%brain_dir = '/Volumes/Expansion Drive/User/Research/Brain_data/BRATS-2/Image_Data/';
%result_path_brain = 'N:\User\Research\Results\textfiles_oct\';
%mask_dir = '/Users/uoft/Dropbox/PhD/BratsAnalysis/CHALLENGE2013_KNN_mha/';
mask_dir = 'N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\';
result_path = 'N:\User\Research\Results\BRATS-2_KNN_withspatialinfo_endDEC29test\';
%result_path = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\KNN_stat_result_for_k\';
% result_path = 'J:\User\Research\Results\BRATS-2_KNN_ENDAUGUST_lowresT2\';
% if exist(result_path,'dir')==0
%     mkdir(result_path)
% end
if exist(result_path_brain,'dir')==0
    mkdir(result_path_brain)
end

Nb_Brains = 30;
mask_type_list = dir(mask_dir);
mask_type_list = mask_type_list(3:end);
 
 
 precision_total = zeros(1,30);
 recall_total = zeros(1,30);
 dice_total = zeros(1,30);
 
 for kkk= 3%result_path_brain = 'J:\User\Research\Results\BRATS-2_KNN_ENDAUGUST_FORUPLOAD\';

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
     for brain_counter = 2:length(mask_list)
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
    
         %% load modalities
        [T1C, T2, FLAIR, truth, info,flair_name] = load_modalities_mac(brain_name_path);
        [height width depth] = size(T1C);
        T1C = double(T1C);
        T2 = double(T2);
        FLAIR = double(FLAIR);
%% reduce the resolution of T2 (check for Maxime)
%        T2 = rescalingT2(T2);
        
        %% creating the mask
        %MASK = make_mask(T1C, T2, FLAIR);
        % or load the mask from a path 
         MASK = load_mask(mask_name_path);
        %% make space and selected points space
       [space,mask_idx,truth_idx]  = make_space(T1C, FLAIR, T2,truth, MASK);
        

        %make_txt_file(space,truth_idx,mask_idx,type, name);
         background = find(sum(space(1:3,:))< 0.001); 
         space(:,background) = [];
         mask_idx(background) = [];
        
         [selected_space, mask_idx] = make_selected_space(space, mask_idx, Nb_healthy,Nb_Selected );
% 
         %% KNN search
%selected_space1 = selected_space(4:6,:);
%space1 = space(4:6,:);
%selected_space = selected_space(1:3,:);
%space = space(1:3,:);
%[IDX,D] = KNN_search(selected_space,space, k);
[IDX,D] = knnsearch(selected_space',space','K', kkk ,'NSMethod','kdtree' ,'Distance','euclidean');

%% assign labels
%selected_space = [selected_space ; selected_space1];
%space = [ space;space1];
segmented = assign_labels(IDX,D,mask_idx, space, height, width, depth);

  
         %% apply median filter
        segmented = medfilt3(segmented, [5,5,5]);
         
         %segmented = u_medfilt3(segmented,5);
         
         %% compute statistics 
         [dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented,truth);
 
  
 %precesion_recall(1,INDEX) = precision ;        
  %precesion_recall(2,INDEX) = recall    ;
  INDEX = INDEX+1;
 t_Elapsed = toc(t_start);
 results_path_root = [result_path,type,'\',name];
 if exist(results_path_root,'dir')==0
    mkdir(results_path_root)
 end
  save([results_path_root,'\','stat_for_',type,'_',name,'.mat'],'precision','recall','dice','jaccard');
% % save the mha file for brain name
%mha_result_name = ['segmented_',type,'_',name,'.mha'];
%mha_result_path = [results_path_root,'\',mha_result_name];
segmented=uint16(segmented);
%writemetaimagefile(mha_result_path, segmented, info.PixelDimensions,info.Offset);
%  
% %%save the mask in MHA
% % MASK = uint16(MASK);
% % mha_MASK_path = [results_path_root,'\','MASK_',type,'_',name,'.mha'];
% writemetaimagefile(mha_MASK_path, MASK, info.PixelDimensions,info.Offset);
 mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
 mha_result_path = [result_path_brain,mha_result_name];
 writemetaimagefile(mha_result_path, segmented, info.PixelDimensions,info.Offset);
 disp(['processing time(load, processes, save) for brain ',type,'_',name,' is ',num2str(t_Elapsed)])
%         
make_txt_file_stats(precision,recall,dice,jaccard, result_path,[type,'_',name])
       
     end
 end
 %% save statistical results
 %save([result_path,'stat_for_k_',num2str(kkk),'.mat'],'precesion_recall');
 
 end
