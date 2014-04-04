clear all
W=5;
Nb_Selected = 500;
Nb_healthy = 500;
dir_name ='N:\BRATS-2\Image_Data\';
 root = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\CHALLENGE2013_KNN_mha\';
 result_path = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\KNN_stat_result_for_k\';
 Nb_Brains = 30;
 root1_list = dir(root);
 root1_list = root1_list(3:end);
 
 precision_total = zeros(1,30);
 recall_total = zeros(1,30);
 dice_total = zeros(1,30);
 
 for kkk= 5:5:50
 INDEX = 1;
 precesion_recall = zeros(2,Nb_Brains);
 for subfoldercounter=1:length(root1_list)
     subforlder1 = root1_list(subfoldercounter);
     root2_list = dir([root,subforlder1.name]);
     root2_list = root2_list(3:end);
     truth_root2_list = dir([dir_name,subforlder1.name]);
     truth_root2_list = truth_root2_list(2:end);
     
     for brain_counter = 1:length(root2_list)
         subforlder2 = root2_list(brain_counter);
         root3_list = dir([root,subforlder1.name,'\',subforlder2.name]);
         root3_list = root3_list(3:end);
         truth_root3_list = dir([dir_name,subforlder1.name,'\',subforlder2.name]);
         truth_root3_list = truth_root3_list(3:end);
         disp(['Processing brain ',subforlder1.name,'_',subforlder2.name, ' with k = ', num2str(kkk)])
         t_start = tic;
         for k = 1:length(root3_list)
             filename_mha = root3_list(k);
        
             if ~isempty(strfind(filename_mha.name,'MASK.nii'))
                seg_path =  [root,subforlder1.name,'\',subforlder2.name,'\',filename_mha.name];
                mha_name = filename_mha.name;
                 MASK = load_nii(seg_path, [], 1);
                 MASK=MASK.img; 
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
            if ~isempty(strfind(filename.name,'T2'))
                filepath = [dir_name,subforlder1.name,'\',subforlder2.name,'\',filename.name];
                files = dir(filepath);
                files = files(3:end);
                for kk = 1: length(files)
                    modality=files(kk);
                    if ~isempty(strfind(modality.name,'T2'))
                        modality_path =  [dir_name,subforlder1.name,'\',subforlder2.name,'\',filename.name,'\',modality.name];
                        T2 = mha_read_volume(modality_path);
                        info  = mha_read_header(modality_path);
                    end
                end 
            end  
            if ~isempty(strfind(filename.name,'Flair'))
                filepath = [dir_name,subforlder1.name,'\',subforlder2.name,'\',filename.name];
                files = dir(filepath);
                files = files(3:end);
                for kk = 1: length(files)
                    modality=files(kk);
                    if ~isempty(strfind(modality.name,'Flair'))
                        modality_path =  [dir_name,subforlder1.name,'\',subforlder2.name,'\',filename.name,'\',modality.name];
                        FLAIR = mha_read_volume(modality_path);
                        info  = mha_read_header(modality_path);
                    end
                end 
            end  
            if ~isempty(strfind(filename.name,'T1c'))
                filepath = [dir_name,subforlder1.name,'\',subforlder2.name,'\',filename.name];
                files = dir(filepath);
                files = files(3:end);
                for kk = 1: length(files)
                    modality=files(kk);
                    if ~isempty(strfind(modality.name,'T1c'))
                        modality_path =  [dir_name,subforlder1.name,'\',subforlder2.name,'\',filename.name,'\',modality.name];
                        T1C = mha_read_volume(modality_path);
                        info  = mha_read_header(modality_path);
                    end
                end 
            end  
            
            
            
         end      
                   
[height width depth] = size(T1C);

%     T1C= T1C(:);
%     T2 = T2(:);
T1C = double(T1C);
T2 = double(T2);
FLAIR = double(FLAIR);
T1C = T1C/max(max(max(T1C)));
T2 = T2/max(max(max(T2)));
FLAIR = FLAIR/max(max(max(FLAIR)));
truth(truth<0.5)=0;                 %healthy 
truth(0.5<truth & truth<1.5)=2;     %tomur
truth(1.5<truth & truth<2.5)=1;     %edema




%% creating the 6 D space (T1,T2,FLAIR,x,y,z)

space = zeros(6,length(T1C(:)));
mask_idx = zeros(1,size(space,2));
% [width, height, depth] = size(T1C);
 
for ROW = 1:height
   for COL=1:width
       for DEP = 1:depth
            
                space(:,depth*(ROW* width + COL)+ DEP)=[T1C(ROW,COL,DEP);T2(ROW,COL,DEP);FLAIR(ROW,COL,DEP);ROW/height;COL/width;DEP/depth];
                mask_idx(depth*(ROW* width + COL)+ DEP) = MASK(ROW,COL,DEP);
           
        end
    end
end

background = find(sum(space(1:3,:))< 0.000001);
space(:,background) = [];
mask_idx(background) = [];
%     for i=1:size(space,1)
%         space(i,:) = space(i,:)/max(space(i,:));
%     end
% idx_downsample = downsample(1:size(space,2), downsamplerate);
% space = space(:,idx_downsample);
selected_points = find(mask_idx> 0 & mask_idx<10);
healthy_points = find(mask_idx ==10);
iindex = downsample(selected_points, ceil(length(selected_points)/Nb_Selected));%sort(uint32(rand(Nb_Selected,1)*length(selected_points)));
h_index = downsample(healthy_points, ceil(length(healthy_points)/Nb_healthy));
selected_space = [space(:,iindex),space(:,h_index)];
mask_idx = mask_idx([iindex,h_index]);
%mask_idx = mask_idx(selected_points);


%% KNN 

% space = space(:,:);   
% space_repmat = repmat(space,[1 1 size(selected_space,2)]);
% selected_space_3d = reshape(selected_space,size(selected_space,1),1,size(selected_space,2));
% selected_space_repmat = repmat(selected_space_3d,[1 size(space,2) 1]);
importance_map = zeros (2,size(space,2));
for i = 1:size(space,2)
    
point = space(:,i);
tmp = abs(selected_space-repmat(point,1,size(selected_space,2)));
tmp = sum(tmp,1);
[tmp IX] = sort(tmp);
val = tmp(1:kkk);
class = mode(mask_idx(IX(1:kkk)));
importance_map(1,i) = val(1);
importance_map(2,i) = class ;      %signifiying the voxel is more likely to be edema


end



%        for j = 1:size(selected_space,2)
%            p_selected = selected_space(1:5,j);
%            distance = sqrt(sum((p_selected - point).^2));
%            imp = min(imp,distance);
%        end
%        importance_map(i) = imp;
%    end

importance_nii = zeros(size(T1C));   
importance_nii_idx = zeros(size(T1C));   

count = 1;
for count = 1:length(importance_map)

        ROW = uint8(space(4,count) * height);
        COL = uint8(space(5,count) * width);
        DEP = uint8(space(6,count) * depth);
        importance_nii(ROW,COL,DEP) = importance_map(1,count);
        importance_nii_idx(ROW,COL,DEP) = importance_map(2,count);  
  
end 

importance_nii_idx(importance_nii_idx==10)=0;
segmented_nii = importance_nii_idx;
importance_nii_idx = medfilt3(importance_nii_idx, [5,5,5]);
[dice,jaccard,precision,recall] = compute_eval_metrics(importance_nii_idx,truth);

% mha_result_name = strrep(flair_name , 'XX.O.MR_Flair',['Seg_',type,'_',name]);
% mha_result_path = [result_path_brain,'\',mha_result_name];
% importance_nii_idx=uint16(importance_nii_idx);
% writemetaimagefile(mha_result_path, importance_nii_idx, infot1c.PixelDimensions,infot1c.Offset);
% MASK = uint16(MASK);
% mha_MASK_path = [result_path_brain,'\','MASK_',type,'_',name,'.mha'];
% writemetaimagefile(mha_MASK_path, MASK, infot1c.PixelDimensions,infot1c.Offset);
 precesion_recall(1,INDEX) = precision ;        
 precesion_recall(1,INDEX) = recall    ;
 INDEX = INDEX+1;
 t_Elapsed = toc(t_start);
 disp(['processing time(load, processes, save) for brain ',subforlder1.name,'_',subforlder2.name,' is ',num2str(t_Elapsed)])
     end
 end
 save([result_path,'stat_for_k_',num2str(kkk),'.mat'],'precesion_recall');
 
 
 end