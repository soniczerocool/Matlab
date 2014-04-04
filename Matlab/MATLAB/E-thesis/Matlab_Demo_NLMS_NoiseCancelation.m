close all
clear all
priv_drawrlsdemo
axis off
signal = sin(2*pi*0.055*(0:1000-1)');

figure
plot(0:199,signal(1:200));
grid; axis([0 200 -2 2]);
title('The information bearing signal');

nvar  = 1.0;                  % Noise variance
noise = randn(1000,1)*nvar;   % White noise
figure
plot(0:999,noise);
title('Noise picked up by the secondary microphone');
grid; axis([0 1000 -4 4]);

nfilt  = fir1(31,0.5);             % 31st order Low pass FIR filter
fnoise = filter(nfilt,1,noise);    % Filtering the noise
d  = signal+fnoise;

figure
plot(0:199,d(1:200));
grid; axis([0 200 -4 4]);
title('Desired input to the Adaptive Filter = Signal + Filtered Noise');

M = 32;                    % Filter order
lam = 1;                   % Exponential weighting factor
delta = 0.1;               % Initial input covariance estimate
w0 = zeros(M,1);           % Initial tap weight vector
P0 = (1/delta)*eye(M,M);   % Initial setting for the P matrix
Zi = zeros(M-1,1);         % FIR filter initial states

% ha=adaptfilt.lms(nord,mu);
% Hadapt = adaptfilt.lms(M,lam,P0,w0,Zi);
mu = 1;                % NLMS step size
Offset = 50;           % NLMS offset
Hadapt =adaptfilt.nlms(32,mu,1,Offset);
Hadapt.PersistentMemory = true;
[y,e] = filter(Hadapt,noise,d);
H = abs(freqz(Hadapt,1,64));
H1 = abs(freqz(nfilt,1,64));

wf = linspace(0,1,64);
plot(wf,H,wf,H1);
xlabel('Normalized Frequency  (\times\pi rad/sample)');
ylabel('Magnitude');
legend('Adaptive Filter Response','Required Filter Response');
grid;
axis([0 1 0 2]);

figure
plot(0:499,signal(1:500),0:499,e(1:500)); grid;
axis([0 500 -4 4]);
title('Original information bearing signal and the error signal');
legend('Original Signal','Error Signal');


error=abs(signal-e);
figure,plot(error);
me=mean(error(800:1000))
se=std(error(800:1000))