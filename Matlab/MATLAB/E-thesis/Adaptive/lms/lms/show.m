
%lms
clear all;
close all;
load anc2_v4;
mu=0.1;
nord=10;
X=convm(x,nord);
[A,E,y]= lms(X,d,mu,nord);
figure
plot(x)
Title('Corelated Noise')
figure
plot(y)
Title('Filter Output')
figure
plot(d)
plot(E)
Title('Expected Speech Signal')



