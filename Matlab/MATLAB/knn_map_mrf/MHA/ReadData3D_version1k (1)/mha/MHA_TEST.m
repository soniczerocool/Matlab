clear all

dir_name = 'N:\BRATS-2\Image_Data\HG\0001\VSD.Brain.XX.O.MR_T1';
list = dir(dir_name);
for j=1:length(list)
    name_modality = list(j);
    if ~isempty(strfind(name_modality.name,'.mha'))
        filepath = [dir_name,'\',name_modality.name];
        V = mha_read_volume(filepath);
        info  = mha_read_header(filepath);
    end
end

V(:,:,100) = 0;

imshow(squeeze(V(:,:,100)),[]);
V=uint16(V);
writemetaimagefile('foo.mha', V, info.PixelDimensions,info.Offset)
%writemetaimagefile(filename, img, resolution, offset)