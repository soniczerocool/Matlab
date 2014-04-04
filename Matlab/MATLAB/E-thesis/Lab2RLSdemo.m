
%lab2 adaptive signal processing
close all;
clear all;
fs=1/1000;
t=0:fs:1;
s=zeros(1,length(t));

fc=10:5:20;

for n=1:length(fc)
   
  A=2*sin(2*pi*fc(n)*t);
  s=s+A;
end

 %now add noise to s
 
 v=rand(1,length(s))
% d=awgn(s,20);
  d=s+v
 SNR=20*log10(norm(s)/norm(v))
  figure(2)
 subplot(4,2,1)
 plot(t,s)
% plot(0:2499,s(1:2500))
title('s(n)without noise')
 subplot(4,2,3)
 plot(t,d)
%  plot(0:2499,v(1:2500))
 title('d(n)after adding noise')
 
 %spectrum of the signals
 
 F_S=abs(fftshift(fft(s)));
 F_D=abs(fftshift(fft(d)));
 subplot(4,2,2)
 plot(F_S)
 subplot(4,2,4)
 plot(F_D)
 
 

 
 
%  


% 

n0_delay=1;
h = dfilt.delay(n0_delay);
delayed_d=filter(h,d);
x=delayed_d;



M = 32;                    % Filter order
lam = 1;                   % Exponential weighting factor
delta = 0.1;               % Initial input covariance estimate
w0 = zeros(M,1);           % Initial tap weight vector
P0 = (1/delta)*eye(M,M);   % Initial setting for the P matrix
Zi = zeros(M-1,1);         % FIR filter initial states



% % 
% filter_length =30;
% beta=.07; 
% X=convm(x,filter_length);
% [A,E,y]= nlms(X,d,beta,filter_length);
%  
%  mu=0.001;
% gama=0.00001;
% nord=10;
% Fs=12000;
% X=convm(x,nord);
% [A,E,y]= llms(X,d,mu,gama,nord);


mu=0.001;
nord=10;
% X=convm(x,nord);
% [A,E,y]= lms(X,d,mu,nord);
%  ha=adaptfilt.lms(nord,mu);
%  [y,E]=filter(ha,delayed_d,d);

Hadapt = adaptfilt.rls(M,lam,P0,w0,Zi);

[y,E] = filter(Hadapt,delayed_d,d);

subplot(4,2,5)
plot(t,y(1:length(t)))
% plot(0:2499,y(1:2500))
subplot(4,2,7)
plot(t,E(1:length(t)))
% plot(0:2499,E(1:2500))
title('error  e(t)');
subplot(4,2,6)
plot(abs(fftshift(fft(y))))
subplot(4,2,8)
plot(abs(fftshift(fft(E))))