 %mex -lut demo2.cpp GCoptimization.cpp graph.cpp LinkedBlockList.cpp maxflow.cpp
 close all
 clear all
 mex -lut markov_network_multilable_general3Dgraph.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\GCoptimization.cpp graph.cpp LinkedBlockList.cpp maxflow.cpp
 I = 1:200;
 I = reshape(I, 5,8,5);
 I = double(I);
%  mask = imread('mask2.pgm');
%  mask = double(mask);
i =1:20;
i = repmat(i,length(I(:)),1);
I_repmat = repmat(I(:),1,size(i,2));
matrix = (abs(I_repmat - i))';

lf = markov_network_multilable_general3Dgraph(I,matrix,500);

 

