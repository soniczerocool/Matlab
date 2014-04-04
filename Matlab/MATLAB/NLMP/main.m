% add to path
addpath('C:\Users\havm2701\Dropbox\PhD\Matlab\MATLAB\NIFTI_20130306')
% load modalities
brain_rootpath = 'N:\User\Research\Brain_data\BRATS-1\Images\';
brain_list = dir(brain_rootpath);




for i=1:length(brain_list)
    subject_test = brain_list(i);
    if ~isempty(strfind(subject_test.name,'SimBRATS'))
        continue;
    end
    if ~isempty(strfind(subject_test.name,'BRATS'))
         brain_test_dir = [brain_rootpath,subject_test.name];
         [T1C_test, T2_test, FLAIR_test, truth_test] = Load_modalites(brain_test_dir);        
    end   
    for j=1:length(brain_list)
        subject_train = brain_list(j);
        if ~isempty(strfind(subject_train.name,'SimBRATS'))|| strcmp(subject_train.name,subject_test.name)
            continue;
        end
        if ~isempty(strfind(subject_train.name,'BRATS')) 
             brain_train_dir = [brain_rootpath,subject_train.name];
             [T1C, T2, FLAIR, truth] = Load_modalites(brain_train_dir);
        display (['Loading ',subject_train.name]);
        end
    
    
  [height,width,depth]= size(T1C) ;
    patch_size = 3;
    window_size = 3;
    for row = 1: height
        for col = 1: width
            for dep = 1:depth
                %patch_test = zeros(nb_features,patch_size.^3);
               
                %patch_test(1,:) = T1c_test();
                sum = 0;
                %window_search = T1C(row-window_size:row+window_size,col-window_size:col+window_size , dep-window_size:dep + window_size);
                for wr = -window_size : window_size
                    for wc = -window_size : window_size
                        for wd = -window_size : window_size
                            
                            for pr = -patch_size : patch_size
                                for pc = -patch_size : patch_size
                                    for pd = -patch_size : patch_size
                                        
                                        dist = abs(T1C_test(min(max(row+pr,1),height), min(max(col+pc,1),width), min(max(dep+pd,1),depth)) - T1C(min(max(row+wr+pr,1),height), min(max(col+wc+pc,1),width), min(max(dep+wd+pd,1),depth)));
                                        sum = sum+ dist;
                                        
                                        
                                    end
                                end
                            end
                            
                        end
                    end
                end
            end
        end
    end
    
    pause;
    end
end

