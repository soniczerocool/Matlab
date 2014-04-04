function segmented = assign_labels(IDX,D,mask_idx,space, height,width,depth)

importance = zeros( height,width,depth);   
segmented = zeros( height,width,depth);   

for count = 1: length(IDX(:,1))

        ROW = uint16(space(4,count) * height);
        COL = uint16(space(5,count) * width);
        DEP = uint16(space(6,count) * depth);
       % importance(ROW,COL,DEP) = min(D(count,:));
        if sum(space(1:3,count)<0.00001)
            segmented(ROW,COL,DEP) = 10;
        else
            segmented(ROW,COL,DEP) = mode(mask_idx(IDX(count,:)));  
        end
end 


segmented(segmented==10)=0;
