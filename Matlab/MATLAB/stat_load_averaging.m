clear all
W=5;
Nb_Selected = 500;
Nb_healthy = 500;
%result_path_brain = 'J:\User\Research\Results\BRATS-2_KNN_ENDAUGUST_FORUPLOAD\';
brain_dir ='J:\User\Research\Results\BRATS-2_KNN_ENDAUGUST_lowresT2\';
mask_dir = 'J:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\';
%result_path = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\KNN_stat_result_for_k\';
result_path = 'J:\User\Research\Results\BRATS-2_KNN_ENDAUGUST_lowresT2\';
if exist(result_path,'dir')==0
    mkdir(result_path)
end
% if exist(result_path_brain,'dir')==0
%     mkdir(result_path_brain)
% end


mask_type_list = dir(mask_dir);
mask_type_list = mask_type_list(3:end);
 
 
 low_p = zeros(29,4);
 low_r = zeros(29,4);
 low_d = zeros(29,4);
 
 full_p = zeros(29,4);
 full_r = zeros(29,4);
 full_d = zeros(29,4);
 
 for kkk= 1:1
 INDEX = 1;
 if kkk==0
     kkk=1;
 end
 

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
         if all(name2 == '0006') & all(type == 'LG')
             continue
         end
         brain_name_path = [brain_type_path,'\',name];
         mask_name_path = [mask_type_path,'\',name];
         TF = strcmp(name,name2);
         if TF == 0
             error('Brain name does not match the mask name!')
         end
          Brain_path = brain_name_path;
         file_list = dir(Brain_path);
        file_list = file_list(3:end);

        for k = 1:length(file_list)
                     filename = file_list(k);
                    if ~isempty(strfind(filename.name,'stat'))
                        low_res_path  = [Brain_path,'\',filename.name];
                        low_res = load(low_res_path);
                    

                        low_p(INDEX,:) = low_res.precision;
                        low_r(INDEX,:) = low_res.recall;
                        low_d (INDEX,:) = low_res.dice;

                    end
        end
         
      
         file_list = dir(mask_name_path);
        file_list = file_list(3:end);

        for k = 1:length(file_list)
                     filename = file_list(k);
                    if ~isempty(strfind(filename.name,'stat'))
                       
                        full_res_path = [mask_name_path,'\',filename.name];
                        full_res  = load(full_res_path);
                        full_p(INDEX,:) = full_res.precision;
                        full_r(INDEX,:) = full_res.recall;
                        full_d (INDEX,:) = full_res.dice;

                    end
        end
         
        
         %% start processing for each brain with id [type,'_',name]
%          t_start = tic;
%         disp(['Processing brain ',type,'_',name, ' with k = ', num2str(kkk)])
    
         %% load modalities
%         [T1C, T2, FLAIR, truth, info,flair_name] = load_modalities(brain_name_path);
%         [height width depth] = size(T1C);
%         T1C = double(T1C);
%         T2 = double(T2);
%         FLAIR = double(FLAIR);
%% reduce the resolution of T2 (check for Maxime)
%         T2 = rescalingT2(T2);
        
        %% creating the mask
        %MASK = make_mask(T1C, T2, FLAIR);
        % or load the mask from a path 
%          MASK = load_mask(mask_name_path);
        %% make space and selected points space
%         [space,mask_idx] = make_space(T1C, FLAIR, T2,truth, MASK);
% 
%         [selected_space, mask_idx] = make_selected_space(space, mask_idx, Nb_healthy,Nb_Selected );
% 
%         %% KNN search
%         % matlab defined function
%         [IDX,D] = knnsearch(selected_space',space','K',kkk,'NSMethod','kdtree' ,'Distance','euclidean');
        % the one I made
        %[IDX,D] = KNN_search(selected_space, space, kkk);
        %% assign labels
%         segmented = assign_labels(IDX,D,mask_idx, space, height, width, depth);

 
        %% apply median filter
%         segmented = medfilt3(segmented, [5,5,5]);
        
        %segmented = u_medfilt3(segmented,5);
        
        %% compute statistics 
%         [dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented,truth);

 
%precesion_recall(1,INDEX) = precision ;        
 %precesion_recall(2,INDEX) = recall    ;
 INDEX = INDEX+1;
%  t_Elapsed = toc(t_start);
%  results_path_root = [result_path,type,'\',name];
%  if exist(results_path_root,'dir')==0
%     mkdir(results_path_root)
%  end
%  save([results_path_root,'\','stat_for_LowResT2_',type,'_',name,'.mat'],'precision','recall','dice','jaccard');
% % save the mha file for brain name
% mha_result_name = ['segmented_LowResT2_',type,'_',name,'.mha'];
% mha_result_path = [results_path_root,'\',mha_result_name];
% segmented=uint16(segmented);
% writemetaimagefile(mha_result_path, segmented, info.PixelDimensions,info.Offset);
 
%%save the mask in MHA
% MASK = uint16(MASK);
% mha_MASK_path = [results_path_root,'\','MASK_',type,'_',name,'.mha'];
% writemetaimagefile(mha_MASK_path, MASK, info.PixelDimensions,info.Offset);
% mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
% mha_result_path = [result_path_brain,mha_result_name];
% writemetaimagefile(mha_result_path, segmented, info.PixelDimensions,info.Offset);
 disp(['processing time(load, processes, save) for brain ',type,'_',name,' is done'])
        

       
     end
 end
 %% save statistical results
 %save([result_path,'stat_for_k_',num2str(kkk),'.mat'],'precesion_recall');
 
 end

 % write a text file to save results
%  Precision = full_p;
%  Recall = full_r;
%  Dice    = full_d;
%  s = struct('Precision',full_p,'Recall',full_r,'Dice', full_d)
  res_path = ['J:\User\Research\Results\T2_full_res_vs_half_res\']
%  dlmwrite(full_res_path,s)
save([res_path , 'T2_full_res_vs_half_res.mat'],'full_r','full_p','low_r','low_p','full_d','low_d' )
full_p_average = sum(full_p)

