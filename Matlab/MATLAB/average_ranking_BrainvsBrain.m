% in this average ranking method, for every brain we find the winning
% method. this is done by getting the rank of each method on every
% statistical measure and see which brain has the lowest. at the end the
% method who has lower rank on all brain wins.
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
Number_methods = 6;
all_metrics= zeros(4,4,3,Number_methods);

for method = 1:Number_methods
    evaluation_matrix_savedir = [all_method_paths{method},'evaluation_matrix\'];
for brain_counter = 1:length(brain_list)
    name = brain_list(brain_counter).name;
    brain_name = [evaluation_matrix_savedir,name];
    stat=load(brain_name);
    all_metrics(:,:,brain_counter,method) = stat.evaluation_matrix;
end
end


nb_brains = 30;
method_average = zeros(nb_brains,Number_methods);
for brain = 1:nb_brains
        per_brain_matrix = all_metrics(:,:,brain,:);
        per_brain_matrix = reshape(per_brain_matrix, 4,4,Number_methods);
        [dummy, sort_matrix] = sort(per_brain_matrix,3,'descend');
        for method = 1:Number_methods    
            sum=0; beta = 0;
            for i=1:4
                for j=1:4
                    
                    value = sort_matrix(i,j,method);
                    if isnan(value)
                        beta = beta+1;
                        continue
                    end
                    sum = sum+value;  
                end
            end
            method_average(brain,method)= sum/(16-beta);
        end

end

[dummy,sort_id] = sort(method_average,2);

[dummy,idx_rank]=sort(mean(sort_id));

method_rank = zeros(1,Number_methods);
 for l=1:Number_methods
     method_rank(idx_rank(l))=l;
 end
% % 

average_f_measure = zeros(nb_methods);
for method_id = 1:nb_methods
    f_measure = zeros(nb_brains);
    for brain_counter = 1:nb_brains
        stat_measure = all_metrics(:,:,brain_counter,method_id);
        Pr = stat_measure(:,1);
        Re = stat_measure(:,2);
        f_measure(brain_counter) = mean(2*(Pr.*Re)./(Pr+Re));
    end
    average_f_measure(method) = sum(f_measure)/nb_brains;
end

save(main_result_dir,method_rank , average_f_measure);