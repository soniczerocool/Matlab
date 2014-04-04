clear all

Nb_Selected = 100;
downsamplerate = 20;
name = 'BRATS_HG0014';
dir_name =['J:\User\Research\Brain_data\BRATS-1\Images\',name];

results_path_root = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\Plots\';
result_path_brain = [results_path_root, name];
    %dir_name =['C:\Users\havm2701\Documents\MATLAB\',name];

list_modality = dir (dir_name);
list_modality = list_modality(2:end);
for j=1:length(list_modality)
    name_modality = list_modality(j);
    if ~isempty(strfind(name_modality.name,'T1C.nii'))
        filepath = [dir_name,'\',name_modality.name];
        T1C = load_nii(filepath, [], 1);
        T1C=T1C.img; 
    end
    if ~isempty(strfind(name_modality.name,'FLAIR.nii'))
        filepath = [dir_name,'\',name_modality.name];
        FLAIR = load_nii(filepath, [], 1);
        FLAIR=FLAIR.img; 
    end
    if ~isempty(strfind(name_modality.name,'T2.nii'))
        filepath = [dir_name,'\',name_modality.name];
        T2 = load_nii(filepath, [], 1);
        T2=T2.img; 
    end
    if ~isempty(strfind(name_modality.name,'truth.nii'))
        filepath = [dir_name,'\',name_modality.name];
        truth = load_nii(filepath, [], 1);
        truth=truth.img; 
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

%% creating the mask


MASK = zeros(size(T1C));
%MASK_slice = imread('Mask_HG14_slice80.png');
decide = false;
while decide == false
    fig2=figure();
    slice = input('enter number of slice (from image J):');
    imshow(T1C(:,:,slice),[])
    MASK_slice = roipoly;
    imshow(MASK_slice);
    MASK(:,:, slice) = MASK_slice;
    d=input('press "1" to continue selecting ROI');
    if d~=1
        break
    end
end

%     imshow
%     MASK_slice(MASK_slice>0)=1;

%% creating the 5 D space (T1,T2,FLAIR,x,y,z)
count = 1;
space = zeros(6,length(T1C(:)));
mask_idx = zeros(1,size(space,2));
% [width, height, depth] = size(T1C);
 for DEP = 1:depth
    for ROW = 1:height
        for COL=1:width  
        space(:,count)=[T1C(ROW,COL,DEP);T2(ROW,COL,DEP);FLAIR(ROW,COL,DEP);ROW/height;COL/width;DEP/depth];
        mask_idx(count) = MASK(ROW,COL,DEP);
        count = count+1 ;
        end
    end
end

%     for i=1:size(space,1)
%         space(i,:) = space(i,:)/max(space(i,:));
%     end
% idx_downsample = downsample(1:size(space,2), downsamplerate);
% space = space(:,idx_downsample);
importance_map = zeros (1,size(space,2));
selected_points = find(mask_idx>0); 
iindex = downsample(selected_points, floor(length(selected_points)/Nb_Selected));%sort(uint32(rand(Nb_Selected,1)*length(selected_points)));

%iindex = sort(uint32(rand(Nb_Selected,1)*length(selected_points))); %select some random points from the selected ROI
%iindex = myunique(iindex);  %remove repeated indexes
selected_space = space(:,iindex); 

% space = space(:,:);   
% space_repmat = repmat(space,[1 1 size(selected_space,2)]);
% selected_space_3d = reshape(selected_space,size(selected_space,1),1,size(selected_space,2));
% selected_space_repmat = repmat(selected_space_3d,[1 size(space,2) 1]);

for i = 1:size(space,2)
point = space(:,i);
tmp = abs(selected_space-repmat(point,1,size(selected_space,2)));
tmp = sum(tmp,1);
val = min(tmp);
importance_map(i) = val;
end

%        for j = 1:size(selected_space,2)
%            p_selected = selected_space(1:5,j);
%            distance = sqrt(sum((p_selected - point).^2));
%            imp = min(imp,distance);
%        end
%        importance_map(i) = imp;
%    end

importance_nii = zeros(size(T1C));   
count = 1;
for DEP = 1:depth
    for ROW = 1:height
        for COL=1:width

        importance_nii(ROW,COL,DEP) = importance_map(count);
        count = count+1 ;
        end
    end
end   
nii = make_nii(importance_nii);    
save_nii(nii, [result_path_brain,'\','T1cT2importancemap_flair2.nii.gz']);

t = input('Set threshold: ');
importance_nii(importance_nii > t) =1;
importance_nii(importance_nii<=t) = 0;
nii_segmented = make_nii(importance_nii);    
save_nii(nii_segmented, [result_path_brain,'\','Tumormap.nii.gz']);


