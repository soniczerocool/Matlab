function [dice_average,precision_average,recall_average] = eval_metrics_all(segment_path,label_path)


 segment_path = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\CHALLENGE2013\';
% 
 label_path='J:\User\Research\Brain_data\BRATS-1\Images\';


list=dir(segment_path);
list = list(3:end);

precision_total= Inf *ones(1,length(list));
recall_total= Inf *ones(1,length(list));
dice_total = Inf *ones(1,length(list));

for i=1:length(list)
    dir_name =[label_path,list(i).name];
    result_path_brain = [segment_path, list(i).name];
    if exist(result_path_brain,'dir')==0
        mkdir(result_path_brain)
    end

    list_modality = dir (dir_name);
    list_modality = list_modality(2:end);
    
    for j=1:length(list_modality)
        name_modality = list_modality(j);
        if ~isempty(strfind(name_modality.name,'truth.nii'))
            filepath = [dir_name,'\',name_modality.name];
            truth = load_nii(filepath, [], 1);
            truth=truth.img; 
        end 
    end
    
    seg_dir_name =[segment_path , list(i).name];
    
    seg_list_files = dir (seg_dir_name);
    seg_list_files = seg_list_files(2:end);
    
    for j=1:length(seg_list_files)
        name_file = seg_list_files(j);
        if ~isempty(strfind(name_file.name,'segmented_tumor.nii'))
            filepath = [seg_dir_name,'\',name_file.name];
            segmented = load_nii(filepath, [], 1);
            segmented=segmented.img; 
        end 
    end
    
    [dice,jaccard,precision,recall] = compute_eval_metrics(segmented,truth);
    

    precision_total(i) = precision;
    recall_total(i) = recall;
    dice_total(i) = dice;
end

precision_average = mean(precision_total);
recall_average = mean(recall_total);
dice_average = mean(dice_total);