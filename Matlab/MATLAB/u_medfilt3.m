function segmented = u_medfilt3(segmented,w)

[height,width,depth] = size(segmented);

for z = 1:depth
             segmented(:,:,z) = medfilt2( segmented(:,:,z),[w,w]);
end