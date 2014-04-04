 %mex -lut demo2.cpp GCoptimization.cpp graph.cpp LinkedBlockList.cpp maxflow.cpp
 close all
 clear all
 mex -lut boykov_jolly_3d.cpp GCoptimization.cpp graph.cpp LinkedBlockList.cpp maxflow.cpp
 I = imread('boy.pgm');
 
 I = double(I);
 [sy,sx] = size(I);
 image = zeros(size(I,1),size(I,2),2);
 image(:,:,1) = I;
 image(:,:,2) = I;
 mask1 = imread('mask3.pgm');
 mask1 = double(mask1);
 mask = 255*ones(size(image));
 mask(:,:,1) = mask1;
%  mask(:,:,2) = mask1;
%matrix = zeros(length(I(:)),50);
%i =1:256;
%i = repmat(i,length(I(:)),1);
%I_repmat = repmat(I(:),1,size(i,2));
%matrix = ((I_repmat - i).^2)';
segmented = zeros(size(I));
lf = boykov_jolly_3d(image,mask,60);
% for y=1:length(sy)
%     for x=1:length(sx)
%  segmented(sy,sx) = lf(y*sx+x);
%     end
% end
lf = reshape(lf,sy,sx,2);
figure,imshow(lf(:,:,1),[])
figure, imshow(I,[])
 