clear all
W=5;
dir_name ='N:\BRATS-2\Image_Data\';
 root = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\CHALLENGE2013_KNN_mha\';
 result_path = 'C:\Users\havm2701\Desktop\UPLOAD_BRATS2013';
 
 root1_list = dir(root);
 root1_list = root1_list(3:end);
 INDEX = 1;
 precision_total = zeros(1,30);
 recall_total = zeros(1,30);
 dice_total = zeros(1,30);
 
 
 
 for i=1:length(root1_list)
     subforlder1 = root1_list(i);
     root2_list = dir([root,subforlder1.name]);
     root2_list = root2_list(3:end);
     truth_root2_list = dir([dir_name,subforlder1.name]);
     truth_root2_list = truth_root2_list(2:end);
     
     for j = 1:length(root2_list)
         subforlder2 = root2_list(j);
         root3_list = dir([root,subforlder1.name,'\',subforlder2.name]);
         root3_list = root3_list(3:end);
         truth_root3_list = dir([dir_name,subforlder1.name,'\',subforlder2.name]);
         truth_root3_list = truth_root3_list(3:end);
         
         for k = 1:length(root3_list)
             filename_mha = root3_list(k);
        
             if ~isempty(strfind(filename_mha.name,'VSD.Brain'))
                seg_path =  [root,subforlder1.name,'\',subforlder2.name,'\',filename_mha.name];
                mha_name = filename_mha.name;
                seg = mha_read_volume(seg_path);
                info  = mha_read_header(seg_path);
             end    
         end
                   
         for k = 1:length(truth_root3_list)
             filename = truth_root3_list(k);
        
            if ~isempty(strfind(filename.name,'OT'))
                filepath = [dir_name,subforlder1.name,'\',subforlder2.name,'\',filename.name];
                files = dir(filepath);
                files = files(3:end);
                for kk = 1: length(files)
                    modality=files(kk);
                    if ~isempty(strfind(modality.name,'OT'))
                        modality_path =  [dir_name,subforlder1.name,'\',subforlder2.name,'\',filename.name,'\',modality.name];
                        truth = mha_read_volume(modality_path);
                        info  = mha_read_header(modality_path);
                    end
                end
            end  
         end      
                   
                  seg = medfilt3(seg,[W,W,W]);
                [dice,jaccard,precision,recall] = compute_eval_metrics(seg,truth);
                precision_total(INDEX) = precision;
                recall_total(INDEX) = recall;
                dice_total(INDEX) = dice;
                seg = uint16(seg);
                medfilt_seg_path = [root,subforlder1.name,'\',subforlder2.name,'\Seg_medfilt3_',subforlder1.name,'_',filename_mha.name];
                writemetaimagefile(medfilt_seg_path, seg, info.PixelDimensions,info.Offset);
                upload_path =[result_path,'\',mha_name];
                save([root,subforlder1.name,'\',subforlder2.name,'\','medfiltered_seg_stats.mat'],'precision','recall','dice')
                writemetaimagefile(upload_path, seg, info.PixelDimensions,info.Offset);
              INDEX = INDEX +1;
           
    end
 end

 