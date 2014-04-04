clear all
close all;
addpath(genpath('C:\Users\havm2701\Dropbox\PhD\Matlab\MATLAB\MHA'))
run('C:\Users\havm2701\Dropbox\PhD\Matlab\MATLAB\vlfeat-0.9.17\toolbox\vl_setup');
global T1 T1C T2 FLAIR 

z_size = zeros(1,30)
compute_statistics = 1;


%
%  brain_name_path = ['N:\User\Research\Brain_data\BRATS-2\Image_Data\',type,'\',name];
%  mask_name_path = ['N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\',type,'\',name];

if compute_statistics
    brain_dir ='N:\User\Research\Brain_data\BRATS-2\Image_Data\';
    % mask_dir = ['N:\User\Research\Brain_data\BRATS-2_MASKS\Image_Data\'];
elseif compute_statistics==0
    brain_dir = 'N:\User\Research\Brain_data\BRATS2013_CHALLENGE\Challenge\';
    mask_dir = 'N:\User\Research\Brain_data\BRATS2013_CHALLENGE_MASKS\Image_Data\';
end
    

 

Nb_Brains = 30;
brain_type_list = dir(brain_dir);
brain_type_list = brain_type_list(3:end);
 for subfoldercounter=1:length(brain_type_list)
     type = brain_type_list(subfoldercounter).name;
     brain_type_path = [brain_dir,type];
     brain_list = dir(brain_type_path);
     brain_list =brain_list(3:end);
   
     for brain_counter =1:length(brain_list)
 
         
         name = brain_list(brain_counter).name;
         brain_name_path = [brain_type_path,'\',name];


t1 = tic();
%% load modalities
        [T1C, T2, FLAIR, truth, T1, info,flair_name] = load_modalities(brain_name_path,compute_statistics);
        [height width depth] = size(T1C);
         T1C = uint16(T1C);
        T2 = uint16(T2);
        FLAIR = uint16(FLAIR);
        T1 = uint16(T1);
        % take backup of the modalities
         truth = uint8(truth);

    

%% Detect the borders of the brain

[ymin,ymax,xmin,xmax,zmin,zmax] = find_boundingbox_borders(T1C);



%% 
% 
disp('***********************************')
disp(['Processing brain ',type,'_',name])
disp('***********************************')

%% write in the textfile





%%
T1C = T1C(xmin:xmax , ymin:ymax , zmin:zmax);
T2 = T2(xmin:xmax , ymin:ymax , zmin:zmax);
FLAIR = FLAIR(xmin:xmax , ymin:ymax , zmin:zmax);
T1 = T1(xmin:xmax , ymin:ymax , zmin:zmax);

if compute_statistics
    truth = truth(xmin:xmax , ymin:ymax , zmin:zmax);
end 

% T1C = T1C/max(max(max(T1C)));
% T2 = T2/max(max(max(T2)));
% FLAIR = FLAIR/max(max(max(FLAIR)));
% truth = truth/max(max(max(truth)));


%% saving the images
brain_name = [type,'_',name];
main_result_root = 'N:\User\Research\Results\ConvNet\images\';
for slice = 1 : size(T1C,3)
    image_dir = [main_result_root,brain_name,'\slice_',num2str(slice),'\'];
    if exist(image_dir,'dir')==0
         mkdir(image_dir) 
    end
    image_t1 = T1(:,:,slice);
    image_t1c = T1C(:,:,slice);
    image_t2 = T2(:,:,slice);
    image_flair = FLAIR(:,:,slice);
    image_truth = truth(:,:,slice);
    imwrite(image_t1,[image_dir,'\T1.pgm'])
    imwrite(image_t1c,[image_dir,'\T1c.pgm'])
    imwrite(image_t2, [image_dir,'\T2.pgm'])
    imwrite(image_flair, [image_dir,'\FLAIR.pgm'])
    imwrite(image_truth, [image_dir,'\truth.pgm'])
I=zeros(size(image_flair,1),size(T1C,2),3);
I = uint16(I);
I(:,:,1)=image_t1c;
I(:,:,2)=image_t2;
I(:,:,3)=image_flair;

imwrite(I,[image_dir,'\RGB_mri.ppm'])    

zip([main_result_root,brain_name,'\slice_',num2str(slice),'.zip'],image_dir)
image_dir(end) = [];
rmdir(image_dir,'s')

end



     end
 end