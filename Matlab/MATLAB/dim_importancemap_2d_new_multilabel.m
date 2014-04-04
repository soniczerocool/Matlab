clear all

Nb_Selected = 500;
downsamplerate = 20;
name = 'BRATS_HG0001';
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
truth(truth<0.5)=0;                 %healthy 
truth(0.5<truth & truth<1.5)=1;     %tomur
truth(1.5<truth & truth<2.5)=2;     %edema


%% creating the mask


MASK = zeros(size(T1C));
%MASK_slice = imread('Mask_HG14_slice80.png');
decide = false;
while decide == false
    fig2=figure();
    sliceT = input('enter slice number for tumor (from image J):');
    
    imshow(T1C(:,:,sliceT),[])
    MASK_sliceT = roipoly;
    MASK_sliceT = 2 * MASK_sliceT;
    imshow(MASK_sliceT);
   
    sliceE = input('enter slice number for edema:');
    imshow(T1C(:,:,sliceE),[])
    MASK_sliceE = roipoly;
    imshow(MASK_sliceE);
    if sliceT == sliceE
        MASK(:,:, sliceT)= min((MASK_sliceT +MASK_sliceE),2);
    else
        MASK(:,:, sliceT) = MASK_sliceT;
        MASK(:,:, sliceE) = MASK_sliceE; 
    end
    
    
    d=input('press "1" to continue selecting ROI');
    if d~=1
        break
    end
end

%     imshow
%     MASK_slice(MASK_slice>0)=1;

%% creating the 6 D space (T1,T2,FLAIR,x,y,z)
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
mask_nii = make_nii(MASK);
save_nii(mask_nii, [result_path_brain,'\','mask.nii.gz']);
selected_points = find(mask_idx>0); 
iindex = downsample(selected_points, floor(length(selected_points)/Nb_Selected));%sort(uint32(rand(Nb_Selected,1)*length(selected_points)));
selected_space = space(:,iindex);
mask_idx = mask_idx(iindex);

edema_points = find(mask_idx==1);
tumor_points = find(mask_idx==2);
%iindex = sort(uint32(rand(Nb_Selected,1)*length(selected_points))); %select some random points from the selected ROI
%iindex = myunique(iindex);  %remove repeated indexes
selected_spaceT = selected_space(:,tumor_points);
selected_spaceE = selected_space(:,edema_points);

% space = space(:,:);   
% space_repmat = repmat(space,[1 1 size(selected_space,2)]);
% selected_space_3d = reshape(selected_space,size(selected_space,1),1,size(selected_space,2));
% selected_space_repmat = repmat(selected_space_3d,[1 size(space,2) 1]);
importance_map = zeros (2,size(space,2));
for i = 1:size(space,2)
point = space(:,i);
tmpT = abs(selected_spaceT-repmat(point,1,size(selected_spaceT,2)));
tmpT = sum(tmpT,1);
valT = min(tmpT);

tmpE = abs(selected_spaceE-repmat(point,1,size(selected_spaceE,2)));
tmpE = sum(tmpE,1);
valE = min(tmpE);

if valE < valT
    importance_map(1,i) = valE;
    importance_map(2,i) = 1;      %signifiying the voxel is more likely to be edema
elseif valT <= valT
    importance_map(1,i) = valT;
    importance_map(2,i) = 2 ;     %signifiying the voxel is more likely to be tumor
end

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
for DEP = 1:depth
    for ROW = 1:height
        for COL=1:width

        importance_nii(ROW,COL,DEP) = importance_map(1,count);
        importance_nii_idx(ROW,COL,DEP) = importance_map(2,count);  
        count = count+1 ;
        end
    end
end   
nii = make_nii(importance_nii);    
save_nii(nii, [result_path_brain,'\','T1cT2importancemap_flair2.nii.gz']);
segmented_nii = importance_nii; 
t = input('Set threshold: ');
segmented_nii(segmented_nii > t) =1;
segmented_nii(segmented_nii<=t) = 0;
segmented_nii = ceil(-(segmented_nii - 0.5));
idx = find(segmented_nii>0);
segmented_nii(idx) = importance_nii_idx(idx);
[dice,jaccard,precision,recall] = compute_eval_metrics(segmented_nii,truth);
segmented_nii = make_nii(segmented_nii);    
save_nii(segmented_nii, [result_path_brain,'\','segmented_tumor.nii.gz']);