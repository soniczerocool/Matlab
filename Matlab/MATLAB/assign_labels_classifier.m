function segmented = assign_labels_classifier(output,mask_idx,space, height,width,depth)

  
segmented = zeros( height,width,depth);   

for count = 1: length(output)

        ROW = uint16(space(4,count) * height);
        COL = uint16(space(5,count) * width);
        DEP = uint16(space(6,count) * depth);
       % importance(ROW,COL,DEP) = min(D(count,:));
        if sum(space(1:3,count)<0.00001)
            segmented(ROW,COL,DEP) = 0;
        else
            segmented(ROW,COL,DEP) = output(count);  
        end
end 


