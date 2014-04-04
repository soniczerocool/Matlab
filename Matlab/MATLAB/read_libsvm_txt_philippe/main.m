close all
clear all


%input_dir = 'N:\User\Research\Results\SVM_OUTPUTS_fromPhilippe\2013_12_17_libsvm_results\libsvm_results_2014\';
input_dir = 'N:\User\Research\Results\MICCAI_2014\SVM\Probabilities\';
result_save_dir = 'N:\User\Research\Results\MICCAI_2014\SVM\Segmented\';
if exist(result_save_dir,'dir')==0
    mkdir(result_save_dir)
end

list_file = dir(input_dir);
list_file = list_file(3:end);

for count = 1:length(list_file)
    name = list_file(count).name;
    inputfile = [input_dir,name];   
    str_list = strread(name,'%s','delimiter', '_');
    type = str_list{1};
    name = str_list{2};
   brain_name_path = ['N:\User\Research\Brain_data\BRATS-2\Image_Data\',type,'\',name];
   [T1C, T2, FLAIR, truth, info,flair_name] = load_modalities_mac(brain_name_path);
   
    mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
    outputfile = [result_save_dir,mha_result_name];
   
    
    posterior_probability = visualize_svm_outputs(inputfile,outputfile,info.Dimensions,info,1,4);
end
