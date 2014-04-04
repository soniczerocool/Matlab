clear all
classifier_name = 'Tree';
result_sub_directory = 'Newmasks_T1c_T2_flair';
main_result_dir = ['N:\User\Research\Results\MICCAI_2014\',classifier_name,'\',result_sub_directory,'\'];
result_save_dir = 'N:\User\Research\Results\BRATS-2_ICPRoldmasks\';

mrf_result_save_dir =  [main_result_dir,'mrf\nodenoising\'];
mrf_denoised_result_save_dir =  [main_result_dir,'mrf\Denoised\'];

boykov_result_save_dir =  [main_result_dir,'boykov\nodenoising\'];
boykov_denoised_result_save_dir =  [main_result_dir,'boykov\Denoised\'];

bykov_ndt_result_save_dir =  [main_result_dir,'boykov_nodataterm\nodenoising\'];
bykov_ndt_denoised_result_save_dir =  [main_result_dir,'boykov_nodataterm\Denoised\'];


knn_medianfilter_result_save_dir =[main_result_dir,'knn_medianfilter\'];
knn_result_save_dir =[main_result_dir,'knn\'];

evaluation_matrix_savedir = [knn_result_save_dir,'evaluation_matrix\'];

brain_list = dir(evaluation_matrix_savedir);
brain_list = brain_list(3:end);

all_method_paths = {knn_result_save_dir;
                    knn_medianfilter_result_save_dir;
                    mrf_denoised_result_save_dir;
                    mrf_result_save_dir;                
                    boykov_denoised_result_save_dir;
                    boykov_result_save_dir;
                    bykov_ndt_result_save_dir;
                    bykov_ndt_denoised_result_save_dir;
                    };
%% 
Number_methods = 5;
nb_brains = 3
all_metrics= zeros(4,4,30,Number_methods);

for method = 1:Number_methods
    evaluation_matrix_savedir = [all_method_paths{method},'evaluation_matrix\'];
for brain_counter = 1:length(brain_list)
    name = brain_list(brain_counter).name;
    brain_name = [evaluation_matrix_savedir,name];
    stat=load(brain_name);
    all_metrics(:,:,brain_counter,method) = stat.evaluation_matrix;
end
end
%avg = mean(all_metrics,3);
method_average = zeros(4,4,Number_methods);
for method=1:Number_methods
    for i=1:4
        for j=1:4
            sum=0; beta = 0;
            for brain = 1:nb_brains
                value = all_metrics(i,j,brain,method);
                if isnan(value)
                    beta = beta+1;
                    continue
                end
                sum = sum+value;  
            end
            method_average(i,j,method)= sum/(nb_brains-beta);
        end
    end
end


[dummy,sort_matrix]=sort(method_average,3,'descend');
ranking_matrix = zeros(4,4,Number_methods);
for i=1:4
   for j=1:4
       rank = sort_matrix(i,j,:);
       for method = 1:Number_methods
        ranking_matrix(i,j,rank(method)) = method;   
           
       end
   end
end

average_rank = zeros(1,Number_methods);

for method = 1:Number_methods
    average_rank(method) = mean2(ranking_matrix(:,:,method));
end
[dummy,idx_rank]= sort(average_rank);
 method_rank_1 = zeros(1,Number_methods);

 for l=1:Number_methods
     method_rank_1(idx_rank(l))=l
 end
