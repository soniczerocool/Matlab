mainfolder='J:\User\Research\Brain_data\BRATS-1\Images\';
results_path_root = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\Plots\';
if exist(results_path_root,'dir')==0
    mkdir(results_path_root)
end

list=dir(mainfolder);
list = list(3:end);
i = 1; 
    dir_name =['J:\User\Research\Brain_data\BRATS-1\Images\',list(i).name];
    result_path_brain = [results_path_root, list(i).name];
    if exist(result_path_brain,'dir')==0
        mkdir(result_path_brain)
    end
    list_modality = dir (dir_name);
    list_modality = list_modality(2:end);
    for j=1:length(list_modality)
        name_modality = list_modality(j);
        if ~isempty(strfind(name_modality.name,'T1C.nii'))
            filepath = [dir_name,'\',name_modality.name];
            T1C = load_nii(filepath, [], 1);
            T1C=T1C.img; 
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
    count = 1;
    space = zeros(5,length(T1C(:)));
    [width, height, depth] = size(T1C);
    for COL=1:width
        for ROW = 1:height
            for DEP = 1:depth
                space(:,count)=[T1C(COL,ROW,DEP);T2(COL,ROW,DEP);COL;ROW;DEP];
                count = count+1 ;
            end
        end
    end
    
    
    
    
    