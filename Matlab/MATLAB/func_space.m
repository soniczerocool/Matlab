function [space,mask_idx] = make_space(T1C, FLAIR, T2,truth, MASK) 



T1C = double(T1C);
T2 = double(T2);
FLAIR = double(FLAIR);
T1C = T1C/max(max(max(T1C)));
T2 = T2/max(max(max(T2)));
FLAIR = FLAIR/max(max(max(FLAIR)));
truth(truth<0.5)=0;                 %healthy 
truth(0.5<truth & truth<1.5)=2;     %tomur
truth(1.5<truth & truth<2.5)=1;     %edema




%% creating the 6 D space (T1,T2,FLAIR,x,y,z)

space = zeros(6,length(T1C(:)));
mask_idx = zeros(1,size(space,2));
% [width, height, depth] = size(T1C);
 
for ROW = 1:height
   for COL=1:width
       for DEP = 1:depth
            
                space(:,depth*(ROW* width + COL)+ DEP)=[T1C(ROW,COL,DEP);T2(ROW,COL,DEP);FLAIR(ROW,COL,DEP);ROW/height;COL/width;DEP/depth];
                mask_idx(depth*(ROW* width + COL)+ DEP) = MASK(ROW,COL,DEP);
           
        end
    end
end



background = find(sum(space(1:3,:))< 0.001);
space(:,background) = [];
mask_idx(background) = [];