close all;
clear all;
clc;
load anc2_v4;
mu=0.6;
gama=0.00001;
nord=10;
Fs=12000;
X=convm(x,nord);
[A,E,y]= llms(X,d,mu,gama,nord);

figure
plot(X)
Title('Corelated Noise')

figure
plot(d)
hold
plot(E,'r')
Title('Expected Speech Signal at mu = 0.6, E')
hold


figure
plot(y)
Title('Filter Coefficient y')


sound(E,Fs)
sound([d(:);E(:)],Fs)
sound([x(:);E(:)],Fs)
