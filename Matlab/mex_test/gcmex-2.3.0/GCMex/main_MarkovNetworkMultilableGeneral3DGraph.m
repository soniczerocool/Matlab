 %mex -lut demo2.cpp GCoptimization.cpp graph.cpp LinkedBlockList.cpp maxflow.cpp
 close all
 clear all
 mex -lut C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\markov_network_multilable5_old.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\GCoptimization.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\graph.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\LinkedBlockList.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\maxflow.cpp
 II = imread('butterflyNoise.pgm');
 [sy,sx] = size(II);
 I1 = double(II);
 I2 = I1;
 I = [I1(:);I2(:)];
 I = reshape(I,sy,sx,2);
%  mask = imread('mask2.pgm');
%  mask = double(mask);
matrix = zeros(length(I(:)),50);
i =1:256;
i = repmat(i,length(I(:)),1);
I_repmat = repmat(I(:),1,size(i,2));
matrix = ((I_repmat - i).^2)';

lf = markov_network_multilable5_old(I,matrix,100);
%imshow(lf,[])
%figure, imshow(I,[])
lf = reshape(lf,sy,sx,2);
imshow(lf(:,:,1),[])