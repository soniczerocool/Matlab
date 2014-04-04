function segmented_volume = label_to_matrix(segmented,space,h,w,d)
height = h;
width = w;
depth = d;
segmented_volume = zeros( height,width,depth);   

for count = 1: length(segmented)

        ROW = uint8(space(4,count) * height);
        COL = uint8(space(5,count) * width);
        DEP = uint8(space(6,count) * depth);
       
        segmented_volume(ROW,COL,DEP) =segmented(count);  
  
end 



