%nlms
clear all;
clc;
close all
load anc2_v4;
beta=1;
nord=10;
Fs=12000;
X=convm(x,nord);
[A,E,y]= nlms(X,d,beta,nord);
figure
plot(X)
Title('Corelated Noise')
figure
plot(y)
Title('Filter Output')
figure
plot(d)
figure
plot(E,'r')
Title('Expected Speech Signal')
hold