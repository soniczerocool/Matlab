 function record_method_evaluation( alpha, beta, save_path3,main_result_dir)
% in this average ranking method, for every brain we find the winning
% method. this is done by getting the rank of each method on every
% statistical measure and see which brain has the lowest. at the end the
% method who has lower rank on all brain wins.

boykov_denoised_result_save_dir =  [main_result_dir,'boykov\Denoised\'];

knn_result_save_dir =[main_result_dir,'knn\'];

evaluation_matrix_savedir = [knn_result_save_dir,'evaluation_matrix\'];

brain_list = dir(evaluation_matrix_savedir);
brain_list = brain_list(3:end);

all_method_paths = {               
                    boykov_denoised_result_save_dir
             
                    };
 
Number_methods = 1;
nb_methods = 1;
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


nb_brains = 30;
method_average = zeros(nb_brains,Number_methods);
for brain = 1:nb_brains
        per_brain_matrix = all_metrics(:,:,brain);
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


average_f_measure = zeros(nb_methods);
for method_id = 1:nb_methods
    sum = 0;beta = 0;
    for brain_counter = 1:nb_brains
        stat_measure = all_metrics(:,:,brain_counter,method_id);
        Pr = stat_measure(:,1);
        Re = stat_measure(:,2);
        for l = 1:4
            value= 2*(Pr(l).*Re(l))./(Pr(l)+Re(l));
           if isnan(value)
                beta = beta+1;
                continue
           end
           sum = sum+value;    
        end
    end
    average_f_measure=sum/4-beta;
end

%save([main_result_dir , 'average_f_measure');
g = fopen(save_path3,'a');
 
    fprintf(g,['1:',num2str(alpha),',2:',num2str(beta),',3:',num2str(average_f_measure)]);
    
    
    fprintf(g,'\n');     
    fclose(g);