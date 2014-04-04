clear all
close all
global T1C T2 FLAIR
Nb_Selected = 500;
k=21;
labels = [0,2,3,4];
Nb_healthy = 500;
downsamplerate = 20;
% name = 'BRATS_HG0001';
% dir_name =['J:\User\Research\Brain_data\BRATS-1\Images\',name];
 results_path_root = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\CHALLENGE2013_KNN_mha\';
  mex -lut markov_network_multilable_Tumor.cpp GCoptimization.cpp graph.cpp LinkedBlockList.cpp maxflow.cpp

 type = 'HG'
 name='0024';
 dir_name = ['J:\User\Research\Brain_data\BRATS-2\Image_Data\',type,'\',name];
 if exist(results_path_root,'dir')==0
    mkdir(results_path_root)
end
result_path_brain = [results_path_root,type,'\', name];
if exist(result_path_brain,'dir')==0
    mkdir(result_path_brain)
end
    %dir_name =['C:\Users\havm2701\Documents\MATLAB\',name];
list_modality_folder = dir (dir_name);
list_modality_folder = list_modality_folder(2:end);
for j=1:length(list_modality_folder)
    name_modality = list_modality_folder(j);
    if ~isempty(strfind(name_modality.name,'T1c'))
        filepath = [dir_name,'\',name_modality.name];
        files = dir(filepath);
        files = files(2:end);
        for kk = 1: length(files)
            modality=files(kk);
            if ~isempty(strfind(modality.name,'T1c'))
                modality_path =  [dir_name,'\',name_modality.name,'\',modality.name]
                T1C = mha_read_volume(modality_path);
                infot1c  = mha_read_header(modality_path);
            end
        end
    end
   
    if ~isempty(strfind(name_modality.name,'T2'))
        filepath = [dir_name,'\',name_modality.name];
        files = dir(filepath);
        files = files(2:end);
        for kk = 1: length(files)
            modality=files(kk);
            if ~isempty(strfind(modality.name,'T2'))
                modality_path =  [dir_name,'\',name_modality.name,'\',modality.name]
                T2 = mha_read_volume(modality_path);
                infot1c  = mha_read_header(modality_path);
            end
        end
    end
    if ~isempty(strfind(name_modality.name,'OT'))
        filepath = [dir_name,'\',name_modality.name];
        files = dir(filepath);
        files = files(2:end);
        for kk = 1: length(files)
            modality=files(kk);
            if ~isempty(strfind(modality.name,'OT'))
                modality_path =  [dir_name,'\',name_modality.name,'\',modality.name]
                truth = mha_read_volume(modality_path);
                infot1c  = mha_read_header(modality_path);
            end
        end
    end
     if ~isempty(strfind(name_modality.name,'Flair'))
        filepath = [dir_name,'\',name_modality.name];
        files = dir(filepath);
        files = files(2:end);
        for kk = 1: length(files)
            modality=files(kk);
            if ~isempty(strfind(modality.name,'Flair'))
                modality_path =  [dir_name,'\',name_modality.name,'\',modality.name];
                flair_name = modality.name;
                FLAIR = mha_read_volume(modality_path);
                infot1c  = mha_read_header(modality_path);
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


slice = input('choose the slice you want to work with ');


%% creating the mask


MASK = zeros(height,width);
while 1
    %MASK_slice = imread('Mask_HG14_slice80.png');
    while 1
        fig2=figure();
         d=input('press "1" to select healthy region');
        if d==1
            

            modality = input('Enter modality to slelect from: ','s');
            imshow_modality(slice,modality);
            MASK_sliceH = roipoly;
            MASK_sliceH = 10 * MASK_sliceH;
            imshow(MASK_sliceH,[]);
            MASK= min((MASK +MASK_sliceH),10);
            
        else 
            break
        end                
    end

    while 1
        d=input('press "1" to select edema region');
        if d==1
            modality = input('Enter modality to slelect from: ','s');
            imshow_modality(slice,modality);
            MASK_sliceE = roipoly;
            MASK_sliceE = 2 * MASK_sliceE;  
            MASK= min((MASK +MASK_sliceE),10);
            imshow(MASK,[]);
        else 
            break
        end 
    end
    
    while 1
         d=input('press "1" to select enhancing tumor');
        if d==1
            
            modality = input('Enter modality to slelect from: ','s');
            imshow_modality(slice,modality);
            MASK_sliceET = roipoly;
            MASK_sliceET = 4* MASK_sliceET;
            MASK= min((MASK +MASK_sliceET),10);
            imshow(MASK,[]);
        else 
            break
        end 
    end

    while 1
        d=input('press "1" to select none_enhancing tumor');
        if d==1
            
            modality = input('Enter modality to slelect from: ','s');
            imshow_modality(slice,modality);
            MASK_sliceNT = roipoly;
            MASK_sliceNT = 3* MASK_sliceNT;
            MASK= min((MASK +MASK_sliceNT),10);
            imshow(MASK);
      else 
            break
        end 
    end


    dd=input('press "1" to continue selecting ROI');
     if dd~=1
          break
     end
end

%     imshow
%     MASK_slice(MASK_slice>0)=1;

%% creating the 5 D space (T1,T2,FLAIR,x,y)
count = 1;
space = zeros(5,length(MASK(:)));
mask_idx = zeros(1,size(space,2));
% [width, height, depth] = size(T1C);
 
 for ROW = 1:height
     for COL = 1:width  
        space(:,count)=[T1C(ROW,COL,slice);T2(ROW,COL,slice);FLAIR(ROW,COL,slice);ROW/height;COL/width];
        mask_idx(count) = MASK(ROW,COL);
        count = count+1 ;
      
    end
end

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
importance_map = zeros (length(labels),size(space,2));
for i = 1:size(space,2)
    
point = space(:,i);
tmp = abs(selected_space-repmat(point,1,size(selected_space,2)));
tmp = sum(tmp,1);
[tmp IX] = sort(tmp);
val = tmp(1:k);
class = mask_idx(IX(1:k));
importance_map(:,i) = hist(class,labels);  % count the number of votes for every class
   
end
black_part = find(sum(space(1:3,:))==0);
importance_map(:,black_part) = 0;
importance_map(1,black_part) = k;
save('importance_map.mat','importance_map');
%        for j = 1:size(selected_space,2)
%            p_selected = selected_space(1:5,j);
%            distance = sqrt(sum((p_selected - point).^2));
%            imp = min(imp,distance);
%        end
%        importance_map(i) = imp;
%    end
%% apply markove network
image = importance_map';
lf = markov_network_multilable_Tumor(image,width,height,21,80);

%%



importance_nii = zeros(size(MASK));   
importance_nii_idx = zeros(size(MASK));   

count = 1;

 for ROW = 1:height
       for COL = 1:width 
        importance_nii_idx(ROW,COL) = lf(count);  
        count = count+1 ;
        end
    end
   
nii = make_nii(importance_nii);
mask_nii = make_nii(MASK);
save_nii(nii, [result_path_brain,'\','6dimportancemap.nii.gz']);
save_nii(mask_nii, [result_path_brain,'\','MASK.nii.gz']);
 
importance_nii_idx(importance_nii_idx==10)=0;
segmented_nii = importance_nii_idx;
importance_nii_idx = medfilt3(importance_nii_idx, [5,5,5]);
[dice,jaccard,precision,recall] = compute_eval_metrics(importance_nii_idx,truth);

segmented_nii = make_nii(importance_nii_idx); 

save_nii(segmented_nii, [result_path_brain,'\',type,'_',name,'_segmented_tumor.nii.gz']);
save([result_path_brain,'\',type,'_',name,'_stat_results.mat'],'precision','recall');

mha_result_name = strrep(flair_name , 'XX.O.MR_Flair',['Seg_',type,'_',name]);
mha_result_path = [result_path_brain,'\',mha_result_name];
importance_nii_idx=uint16(importance_nii_idx);
writemetaimagefile(mha_result_path, importance_nii_idx, infot1c.PixelDimensions,infot1c.Offset);
MASK = uint16(MASK);
mha_MASK_path = [result_path_brain,'\','MASK_',type,'_',name,'.mha'];
writemetaimagefile(mha_MASK_path, MASK, infot1c.PixelDimensions,infot1c.Offset);