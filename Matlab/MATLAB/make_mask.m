function MASK = make_mask(T1C, T2, FLAIR)
display('***** Make New Mask*****')
[height,width,depth] = size(T1C);
MASK = zeros(size(T1C));
while 1
    %MASK_slice = imread('Mask_HG14_slice80.png');
    while 1
        fig2=figure();
         d=input('press "1" to select healthy region');
        if d==1
            sliceH = input('enter slice number for healthy (from image J):');
%             sliceH = depth - sliceH;
            modality = input('Enter modality to slelect from: ','s');
            view = input('Enter view to slelect from: ');
            imshow_modality(sliceH,modality,view);
            
            roi= roipoly;
            if view ~= 3 
                roi = flipud(roi)';
            elseif view ==3
                roi = roi';
            end
            MASK_sliceH = roi;
            MASK_sliceH = 10 * MASK_sliceH;
            imshow(roi);
            roi_index = find(roi>0);

            if view ==1
            MASK_sliceH =MASK(sliceH,:,:);
            MASK_sliceH(roi_index)= 10;
            MASK(sliceH,:,:) = MASK_sliceH;
            elseif view ==2
            MASK_sliceH =MASK(:,sliceH,:);
            MASK_sliceH(roi_index)= 10;
            MASK(:,sliceH,:) = MASK_sliceH;
            elseif view == 3
            MASK_sliceH =MASK(:,:,sliceH);
            MASK_sliceH(roi_index)= 10;
            MASK(:,:,sliceH) = MASK_sliceH;
            end
        else 
            break
        end                
    end
    while 1
        d=input('press "1" to select edema region');
        if d==1
            sliceE = input('enter slice number for edema (from image J):');
%             sliceE = depth - sliceE;
            modality = input('Enter modality to slelect from: ','s');
            view = input('Enter view to slelect from: ');
            imshow_modality(sliceE,modality,view);
             roi= roipoly;
            if view ~= 3 
                roi = flipud(roi)';
            elseif view ==3
                roi = roi';
            end
            imshow(roi);
            roi_index = find(roi>0);

            if view ==1
            MASK_sliceE =MASK(sliceE,:,:);
            MASK_sliceE(roi_index)= 2;
            MASK(sliceE,:,:) = MASK_sliceE;
            elseif view ==2
            MASK_sliceE =MASK(:,sliceE,:);
            MASK_sliceE(roi_index)= 2;
            MASK(:,sliceE,:) = MASK_sliceE;
            elseif view == 3
            MASK_sliceE =MASK(:,:,sliceE);
            MASK_sliceE(roi_index)= 2;
            MASK(:,:,sliceE) = MASK_sliceE;
            end
        else 
            break
        end                
    end

    while 1
         d=input('press "1" to select enhancing tumor');
        if d==1
            sliceET = input('enter slice number for enhancing tumor:');
%             sliceET = depth - sliceET;
            modality = input('Enter modality to slelect from: ','s');
            view = input('Enter view to slelect from: ');
            imshow_modality(sliceET,modality,view);
            roi= roipoly;
            if view ~= 3 
                roi = flipud(roi)';
            elseif view ==3
                roi = roi';
            end
            %MASK_sliceET = roi;
            %MASK_sliceET = 4 * MASK_sliceET;
            imshow(roi);
            roi_index = find(roi>0);

            if view ==1
            MASK_sliceET =MASK(sliceET,:,:);
            MASK_sliceET(roi_index)= 4;
            MASK(sliceET,:,:) = MASK_sliceET;
            elseif view ==2
            MASK_sliceET =MASK(:,sliceET,:);
            MASK_sliceET(roi_index)= 4;
            MASK(:,sliceET,:) = MASK_sliceET;
            elseif view == 3
            MASK_sliceET =MASK(:,:,sliceET);
            MASK_sliceET(roi_index)= 4;
            MASK(:,:,sliceET) = MASK_sliceET;
            end
        else 
            break
        end                
    end

    while 1
        d=input('press "1" to select none_enhancing tumor');
        if d==1
            sliceNT = input('enter slice number for non enhancing tumor:');
%             sliceNT = depth - sliceNT;
            modality = input('Enter modality to slelect from: ','s');
             view = input('Enter view to slelect from: ');

            imshow_modality(sliceNT,modality,view);
            roi= roipoly;
            if view ~= 3 
                roi = flipud(roi)';
            elseif view ==3
                roi = roi';
            end
            MASK_sliceNT = roi;
            MASK_sliceNT = 3 * MASK_sliceNT;
            imshow(roi);
            roi_index = find(roi>0);

            if view ==1
            MASK_sliceNT =MASK(sliceNT,:,:);
            MASK_sliceNT(roi_index)= 3;
            MASK(sliceNT,:,:) = MASK_sliceNT;
            elseif view ==2
            MASK_sliceNT =MASK(:,sliceNT,:);
            MASK_sliceNT(roi_index)= 3;
            MASK(:,sliceNT,:) = MASK_sliceNT;
            elseif view == 3
            MASK_sliceNT =MASK(:,:,sliceNT);
            MASK_sliceNT(roi_index)= 3;
            MASK(:,:,sliceNT) = MASK_sliceNT;
            end
        else 
            break
        end                
    end

    dd=input('press "1" to continue selecting ROI');
     if dd~=1
          break
     end
end