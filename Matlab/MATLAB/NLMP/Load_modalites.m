function [T1C,T2,FLAIR,truth] = Load_modalites(dir_name)
%file_list = dir(brain_dir);

list_modality = dir (dir_name);
list_modality = list_modality(2:end);
for j=1:length(list_modality)
    name_modality = list_modality(j);
    if ~isempty(strfind(name_modality.name,'T1C'))
        filepath = [dir_name,'\',name_modality.name];
        T1C = load_nii(filepath, [], 1);
        T1C=T1C.img; 
    end

    if ~isempty(strfind(name_modality.name,'T2'))
        filepath = [dir_name,'\',name_modality.name];
        T2 = load_nii(filepath, [], 1);
        T2=T2.img; 
    end
    
        if ~isempty(strfind(name_modality.name,'FLAIR'))
        filepath = [dir_name,'\',name_modality.name];
        FLAIR = load_nii(filepath, [], 1);
        FLAIR=FLAIR.img; 
        end

        if ~isempty(strfind(name_modality.name,'truth'))
        filepath = [dir_name,'\',name_modality.name];
        truth = load_nii(filepath, [], 1);
        truth=truth.img; 
        end
end