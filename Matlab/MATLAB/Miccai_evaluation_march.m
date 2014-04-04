clear all
Number_methods = 40;
 iter = 1
path1 = 'N:\User\Research\Results\MICCAI_2014_MARCH\'
dim_dir=dir(path1);
dim_dir=dim_dir(3:end);
s = 1
for i=1:2
    path2 = [path1,dim_dir(i).name];
    classifier_dir = dir(path2);
    classifier_dir = classifier_dir(3:end);
    for j=1:length(classifier_dir)
        path3 = [path2,'\',classifier_dir(j).name,'\'];
        main_result_dir = path3
        
        %%%%%%%%%%%%
        mrf_result_save_dir =  [main_result_dir,'mrf\nodenoising\'];
        mrf_denoised_result_save_dir =  [main_result_dir,'mrf\Denoised\'];

        boykov_result_save_dir =  [main_result_dir,'boykov\nodenoising\'];
        boykov_denoised_result_save_dir =  [main_result_dir,'boykov\Denoised\'];

        %bykov_ndt_result_save_dir =  [main_result_dir,'boykov_nodataterm\nodenoising\'];
       % bykov_ndt_denoised_result_save_dir =  [main_result_dir,'boykov_nodataterm\Denoised\'];


        knn_medianfilter_result_save_dir =[main_result_dir,classifier_dir(j).name,'_medianfilter\'];
        knn_result_save_dir =[main_result_dir,classifier_dir(j).name,'\'];

        evaluation_matrix_savedir = [knn_result_save_dir,'evaluation_matrix\'];

        brain_list = dir(evaluation_matrix_savedir);
        brain_list = brain_list(3:end);

        method_paths = {
                            knn_result_save_dir;
                            knn_medianfilter_result_save_dir;
                            mrf_denoised_result_save_dir;
                            mrf_result_save_dir;                
                            boykov_denoised_result_save_dir;
                            boykov_result_save_dir;
                           % bykov_ndt_result_save_dir;
                          %  bykov_ndt_denoised_result_save_dir;
                            };
        for jj = 1:length(method_paths)
            all_method_paths{s} = method_paths{jj}
            s = s+1
        end
                           
 
       %%%%%%%%%%%%
       iter = iter +1;
    end
end
% 
Number_methods = length(all_method_paths);
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


%  nb_brains = 30;
% method_average = zeros(nb_brains,Number_methods);
% for brain = 1:nb_brains
%         per_brain_matrix = all_metrics(:,:,brain,:);
%         per_brain_matrix = reshape(per_brain_matrix, 4,4,Number_methods);
%         [dummy, sort_matrix] = sort(per_brain_matrix,3,'descend');
%         for method = 1:Number_methods    
%             sum=0; beta = 0;
%             for i=1:4
%                 for j=1:4
%                     
%                     value = sort_matrix(i,j,method);
%                     if isnan(value)
%                         beta = beta+1;
%                         continue
%                     end
%                     sum = sum+value;  
%                 end
%             end
%             method_average(brain,method)= sum/(16-beta);
%         end
% 
% end
% 
% [dummy,sort_id] = sort(method_average,2);
% 
% [dummy,idx_rank]=sort(mean(sort_id));
% 
% method_rank = zeros(1,Number_methods);
%  for l=1:Number_methods
%      method_rank(idx_rank(l))=l;
%  end
% 
% 
nb_brains = 30
average_f_measure = zeros(1,Number_methods);
for method_id = 1:Number_methods
    brains_stat = zeros(1,nb_brains); 
    for brain_counter = 1:nb_brains
        stat_measure = all_metrics(:,:,brain_counter,method_id);
        Pr = stat_measure(:,1);
        Re = stat_measure(:,2);
        sum = 0;beta = 0;
        for l = 1:4
            value= 2*(Pr(l).*Re(l))./(Pr(l)+Re(l));
           if isnan(value)
                beta = beta+1;
                continue
           end
           sum = sum+value;    
        end
        brains_stat(brain_counter)=sum/4-beta;
    end
    average_f_measure(method_id)= mean(brains_stat);
end
