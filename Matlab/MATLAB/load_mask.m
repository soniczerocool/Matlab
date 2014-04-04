function MASK = load_mask(mask_path)
% a function to load the mask form nii or mha file

file_list = dir(mask_path);
file_list = file_list(3:end);
         
         for k = 1:length(file_list)
             filename = file_list(k);
        
             if ~isempty(strfind(filename.name,'MASK'))
                seg_path =  [mask_path,filesep,filename.name];
                [pathstr,name,ext] = fileparts(seg_path)
                switch ext
                    case '.nii.gz'
                        MASK = load_nii(seg_path, [], 1);
                        MASK=MASK.img;
                    case '.mha'
                   MASK =  mha_read_volume(seg_path);
                end
               
             end    
         end