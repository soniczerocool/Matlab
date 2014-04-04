clear all
W=5;
Nb_Selected = 500;
Nb_healthy = 500;
brain_dir ='J:\User\Research\Brain_data\BRATS-2\Image_Data\';
mask_dir = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\CHALLENGE2013_KNN_mha\';
%result_path = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\KNN_stat_result_for_k\';
result_path = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\MASKS';
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
        
        mask_idx = find(MASK>0);
        MASK(mask_idx) = truth(mask_idx);
       mha_MASK_path = [result_path,'\','MASK_',type,'_',name,'.mha'];
writemetaimagefile(mha_MASK_path, MASK, info.PixelDimensions,info.Offset);
     end
 end
 %% save statistical results
 %save([result_path,'stat_for_k_',num2str(kkk),'.mat'],'precesion_recall');
 
 end