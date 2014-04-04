function imshow_modality(sliceE,modality,view)
global T1C
[height,width,depth] = size(T1C);

switch modality
    case 'T1C'
        
        if view ==3
            imshow(T1C(:,:,sliceE)',[])
        elseif view == 1
            imshow(flipud(reshape(T1C(sliceE,:,:),width,depth)'),[])
        elseif view == 2
             imshow(flipud(reshape(T1C(:,sliceE,:),height,depth)'),[])
        end
    case 'T2'
        global T2
        if view == 3
            imshow(T2(:,:,sliceE)',[])
        elseif view == 1
            imshow(flipud(reshape(T2(sliceE,:,:),width,depth)'),[])
        elseif view == 2
            imshow(flipud(reshape(T2(:,sliceE,:),height,depth)'),[])
        end
    case 'FLAIR'
        global FLAIR
        if view ==3
            imshow(FLAIR(:,:,sliceE)',[])
        elseif view == 1
            imshow(flipud(reshape(FLAIR(sliceE,:,:),width,depth)'),[])
        elseif view == 2
             imshow(flipud(reshape(FLAIR(:,sliceE,:),height,depth)'),[])
        end
end

