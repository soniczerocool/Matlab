clear all
close all;
addpath(genpath('C:\Users\havm2701\Dropbox\PhD\Matlab\MATLAB\MHA'))
run('C:\Users\havm2701\Dropbox\PhD\Matlab\MATLAB\vlfeat-0.9.17\toolbox\vl_setup');
mex -lut C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\markov_network_multilable_general3Dgraph.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\GCoptimization.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\graph.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\LinkedBlockList.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\maxflow.cpp
mex -lut C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\interactiveGraphcut_withdataterms.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\GCoptimization.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\graph.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\LinkedBlockList.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\maxflow.cpp
global T1C T2 FLAIR 
%classifier_name = 'Tree';
classifier_name = 'kNN'
% classifier_name = 'SVM';
no_spatial_feature = 0; % set to one if you want to remove the spatial features
k=3;

if ~strcmp(classifier_name,'kNN')
    k = NaN;
end

% output directories
if no_spatial_feature==1
    result_sub_directory = 'Newmasks_T1c_T2_flair'
else
    result_sub_directory = 'Newmasks_6Dfeature'
end


%brain_dir = '/Volumes/Expansion Drive/User/Research/Brain_data/BRATS-2/Image_Data/';
%result_path_brain = 'N:\User\Research\Results\textfiles_oct\';
%mask_dir = '/Users/uoft/Dropbox/PhD/BratsAnalysis/CHALLENGE2013_KNN_mha/';


%result_path = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\KNN_stat_result_for_k\';
% result_path = 'J:\User\Research\Results\BRATS-2_KNN_ENDAUGUST_lowresT2\';

% 
% if exist(textfiles_save_dir,'dir')==0
%     mkdir(textfiles_save_dir)
% end
% 
% save_path1 = [textfiles_save_dir,'summery_evaluation_methods.txt'];
% % if exist(result_path_brain,'dir')==0
% %     mkdir(result_path_brain)
% % end

%% define paths and parameters
Nb_Selected = 1000;

Nb_healthy = 1000;
downsamplerate = 20;




%  brain_name_path = ['N:\User\Research\Brain_data\BRATS-2\Image_Data\',type,'\',name];
%  mask_name_path = ['N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\',type,'\',name];

 mask_dir = ['N:\User\Research\Brain_data\BRATS-2_MASKS\Image_Data\'];
 brain_dir ='N:\User\Research\Brain_data\BRATS-2\Image_Data\';
%mask_dir = 'N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\';
alpha_vector = 0.03:0.01:0.1;
beta_vector = [0.001:0.001:0.01,0.01:0.01:0.1];
for a =1:10; 
    for b =1:20
 

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
%           if type == 'HG'
% %              error('finished')
%                continue;
%           end
t1 = tic();
%% load modalities
        [T1C, T2, FLAIR, truth, info,flair_name] = load_modalities(brain_name_path);
        [height width depth] = size(T1C);
        T1C = double(T1C);
        T2 = double(T2);
        FLAIR = double(FLAIR);
        % take backup of the modalities
        truth_backup = truth;

%% reduce the resolution of T2 (check for Maxime)
%        T2 = rescalingT2(T2);
        
%% creating the mask and ROI
   % MASK = make_mask(T1C, T2, FLAIR);
% or load the mask from a path
%mask_name_path = 'N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\HG\0026';
 %mask_name_path = ['N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\',type,'\',name];
 mask_name_path = [mask_dir,type,'\',name];
 MASK = load_mask(mask_name_path);
 backup_MASK = MASK;
% you can Edit the mask
     %MASK = edit_mask(T1C, T2, FLAIR, MASK);
 
%[ymin,ymax,xmin,xmax,zmin,zmax] = make_regtangle(T1C);  
% save ROI (regtangle) coordinate points
 


%% Detect the borders of the brain

[ymin,ymax,xmin,xmax,zmin,zmax] = find_boundingbox_borders(T1C);


%%% load reg coordinates from the saved location  => this method is
%%% obseleate in the new version. In the new version we detect the borders
%%% of the brain by find_boundingbox_borders function 


% roi_save_dir = ['N:\User\Research\Brain_data\BRATS-2_roi_coordinates\Image_Data\',type,'\',name];
%  roi_save_path = [roi_save_dir,'\ROI_',type,'_',name,'.mat'] ; 
% load(roi_save_path)
%% 
% %  if exist(roi_save_dir,'dir')==0
%     mkdir(roi_save_dir)
%  end
%    save(roi_save_path,'ymin','ymax','xmin','xmax','zmin','zmax');

 
% mask_save_path = [mask_save_dir,'\MASK_',type,'_',name,'.mha'] ;     
% MASK = uint16(MASK);
% writemetaimagefile(mask_save_path, MASK, info.PixelDimensions,info.Offset);     
% MASK = double(MASK);  
disp('***********************************')
disp(['Processing brain ',type,'_',name, ' with k = ', num2str(k)])
disp('***********************************')






%%
T1C = T1C(xmin:xmax , ymin:ymax , zmin:zmax);
T2 = T2(xmin:xmax , ymin:ymax , zmin:zmax);
FLAIR = FLAIR(xmin:xmax , ymin:ymax , zmin:zmax);

MASK = MASK(xmin:xmax , ymin:ymax , zmin:zmax);
truth = truth(xmin:xmax , ymin:ymax , zmin:zmax);
 
%%  make space and selected points space( MADE CHANGES FOR MAKING TEXT FILES >> FIX LATER)
  
[space,mask_idx,truth_idx] = make_space(T1C, FLAIR, T2,truth, MASK);
 space2 = space;
         background = find(sum(space(1:3,:))< 0.001); 
         space(:,background) = [];
         mask_idx(background) = [];
         truth_idx(background) = [];


% Remove the spatial features
if no_spatial_feature ==1
    cordinates = space(4:6,:);
    space = space(1:3,:);
end


[selected_space, mask_idx] = make_selected_space(space, mask_idx, Nb_healthy,Nb_Selected );


%  make_txt_file(space,selected_space,truth_idx,mask_idx,type, name);

  
 [h,w,d] = size(T1C);
% 
%% KNN search classifier
t2 = tic();
kkk=k;
[IDX,D] = knnsearch(selected_space',space','K', kkk ,'NSMethod','kdtree' ,'Distance','euclidean');
knn_time = toc(t2);
%[ind,knndist] = kdtreeKNN(selected_space,mask_idx,space,k);
%IDX = ind';
%D = knndist';

% assign labels
if no_spatial_feature ==1
    space = [space;cordinates];
end
%segmented_volume = assign_labels(IDX,D,mask_idx, space, h, w, d);

% make data-term_KNN

matrix = make_dataterm_matrix(IDX,mask_idx,space,T1C);
matrix = (matrix/k);
%energy = -log(matrix);
energy = 1 - matrix;


% %
%% APPLY graphcut BOYKOV JOLLY with dataterms

        alpha = alpha_vector(a);
        beta = beta_vector(b);
        if beta == 0.004
            continue
        end
        alpha = 0.01
        beta = 0.005
disp('***********************************')
disp(['Processing brain ',type,'_',name, ' with alpha: ', num2str(alpha),'beta: ',num2str(beta)])
disp('***********************************')       
        
 result_dir =['N:\User\Research\Results\MICCAI_2014_find_parameters\',classifier_name,'\',result_sub_directory,'\'];
main_result_dir = [result_dir,'a',num2str(alpha),'b',num2str(beta),'\'];
save_path3 = [result_dir,'parameter_stats.txt'];
textfiles_save_dir =[main_result_dir,'textresult\'];

%alpha = 0.05;
 energy = alpha * energy;
% beta = .0051;
 t3 = tic();
lf = interactiveGraphcut_withdataterms(space2,double(MASK),energy,beta);
graphcut_time = toc(t3)
[h,w,d] = size(T1C);
  lf = reshape(lf,h, w, d);
   segmented_volume = lf ;
   % apply median filter
segmented_volume = medfilt3(segmented_volume, [3,3,3]);
% compute statistics 
segmented_brain = zeros(height,width,depth);
segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;

 segmented_brain(segmented_brain == 3)=4;
 segmented_brain(segmented_brain == 1)=3;
[dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented_brain,truth);
 stats = [precision;recall;dice;jaccard];
avg = sum(stats(:))/16;
%%%%%%%%%%%% save results
% mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
% mha_result_path = [boykov_result_save_dir,mha_result_name];
% segmented_brain = uint16(segmented_brain);
% writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset); 
%%%% save the matfile for statistics
boykov_denoised_result_save_dir =  [main_result_dir,'boykov\Denoised\'];
if exist(boykov_denoised_result_save_dir,'dir')==0
    mkdir(boykov_denoised_result_save_dir)
end
evaluation_matrix_savedir = [boykov_denoised_result_save_dir,'evaluation_matrix\'];
evaluation_matrix_savepath = [evaluation_matrix_savedir,type,'_',name,'.mat'];
evaluation_matrix = [precision;recall;dice;jaccard];
if exist(evaluation_matrix_savedir,'dir')==0
    mkdir(evaluation_matrix_savedir)
end
 save(evaluation_matrix_savepath,'evaluation_matrix','precision','recall','dice','jaccard')   


%  total_time =  toc(t1);
%  display(['Timings for brain: ',type,'_',name,':']) 
%  display(['Total processing: ',num2str(total_time) ])
%  display(['kkn: ',num2str(knn_time)]);
%  display(['graphcut: ',num2str(graphcut_time)])
     end
 end
 % record_method_evaluation(alpha, beta, save_path3, main_result_dir);
 
    end
end