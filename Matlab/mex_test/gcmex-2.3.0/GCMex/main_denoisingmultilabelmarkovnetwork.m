 %mex -lut demo2.cpp GCoptimization.cpp graph.cpp LinkedBlockList.cpp maxflow.cpp
 close all
 clear all
  mex -lut markov_network_multilable_general2Dgraph.cpp GCoptimization.cpp graph.cpp LinkedBlockList.cpp maxflow.cpp
 I = imread('butterflyNoise.pgm');
 I = double(I);
%  mask = imread('mask2.pgm');
%  mask = double(mask);
matrix = zeros(length(I(:)),50);
i =linspace(0,255,50);
i = repmat(i,length(I(:)),1);
I_repmat = repmat(I(:),1,size(i,2));
matrix = (I_repmat - i).^2;
matrix = matrix';
lf = markov_network_multilable_grid2Dgraph(I,100);
%lf = markov_network_multilable_general2Dgraph(I,matrix,500);
imshow(lf,[])
figure, imshow(I,[])
 