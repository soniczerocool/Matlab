function T2_lowres = rescalingT2(T2)

%rescaling T2 for use to verify if lower resolution T2 would give an
%acceptible result in the KNN search method(interactive method).

 

[height width depth] = size(T2);
h = 1:2:height;
w = 1:2:width;
d = 1:2:depth;
T2_downsampled = T2(h,w,d);
T2_lowres = Inf * ones(size(T2));
for i = 1:length(d)
    T2_lowres(:,:,2*i-1) = imresize(T2_downsampled(:,:,i),2,'bilinear');
    T2_lowres(:,:,2*i) = (imresize(T2_downsampled(:,:,i),2,'bilinear') + imresize(T2_downsampled(:,:,min(i+1,length(d))),2,'bilinear'))/2;
end
