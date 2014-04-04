clear all

W=5;
Nb_Selected = 500;
Nb_healthy = 500;
brain_dir ='J:\User\Research\Brain_data\BRATS-2\Image_Data\';
mask_dir = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\CHALLENGE2013_KNN_mha\';
%result_path = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\KNN_stat_result_for_k\';
result_path = 'J:\User\Research\Results\BRATS-2_KNN\';
Nb_Brains = 30;
mask_type_list = dir(mask_dir);
mask_type_list = mask_type_list(3:end);
 
 
 precision_total = zeros(1,30);
 recall_total = zeros(1,30);
 dice_total = zeros(1,30);
 
 for kkk= 1:1
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
     mask_list = mask_list(3:end);   
     for brain_counter = 1:length(mask_list)
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
        [T1C, T2, FLAIR, truth, info] = load_modalities(brain_name_path);
        [height width depth] = size(T1C);
        T1C = double(T1C);
        T2 = double(T2);
        FLAIR = double(FLAIR);

        %% creating the mask
        %MASK = make_mask(T1C, T2, FLAIR);
        % or load the mask from a path 
         MASK = load_mask(mask_name_path);
       
       
        
        %% compute statistics 
        [dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(MASK,truth);

 
%precesion_recall(1,INDEX) = precision ;        
 %precesion_recall(2,INDEX) = recall    ;
 INDEX = INDEX+1;
 t_Elapsed = toc(t_start);
 results_path_root = [result_path,type,'\',name];
 if exist(results_path_root,'dir')==0
    mkdir(results_path_root)
 end
 save([results_path_root,'\','stat_for_mask_',type,'_',name,'.mat'],'precision','recall','dice','jaccard');
% save the mha file for brain name
%mha_result_name = ['segmented',type,'_',name,'.mha'];
%mha_result_path = [results_path_root,'\',mha_result_name];
%segmented=uint16(segmented);
%writemetaimagefile(mha_result_path, segmented, info.PixelDimensions,info.Offset);
 
 disp(['processing time(load, processes, save) for brain ',type,'_',name,' is ',num2str(t_Elapsed)])
        

       
     end
 end
 %% save statistical results
 %save([result_path,'stat_for_k_',num2str(kkk),'.mat'],'precesion_recall');
 
 end
