clear all;
clc;

%=================   an arbitrary signal 
f=50;
L=500;
Ts=1/500;
Fs=1/Ts;
t=Ts*[0:500-1];

signal =2.7*( sin(2*pi*f*t+pi)+(1/3)*sin(7.6*pi*f*t) ) ;

v=1*randn(size(t));
d=signal' + v';

subplot(4,2,1)
plot(0:499,signal(1:500),'g')
title('signal without noise');

subplot(4,2,3)
plot(0:499,d(1:500),'r')
title('signal with noise');


% ===================  Frequency Components

Hs = spectrum.periodogram('Hamming');  %periodogram
NFFT = 2^nextpow2(L); % Next power of 2 from length of signal
AA=psd(Hs,signal,'Fs',Fs,'NFFT',NFFT,'SpectrumType','twosided');
Pow_of_signal = avgpower(AA);

BB=psd(Hs,v,'Fs',Fs,'NFFT',NFFT,'SpectrumType','twosided');
Pow_of_noise = avgpower (BB) ;

SNR = Pow_of_signal/Pow_of_noise;
SNR_DB=10*log(SNR)

subplot(4,2,2)
plot(abs( fftshift(fft(signal) ) ))
title('signal without noise');

subplot(4,2,4)
plot(abs( fftshift(fft(d) ) ))
title('signal with noise (d)');



% ================  to delay the signal
 
n0_delay=3;
h = dfilt.delay(n0_delay);
delayed_d=filter(h,d);

filter_length =5;
miu= 0.006;

% ha=adaptfilt.lms(filter_length,miu);
% [y,e]=filter(ha,delayed_d,d);

x=delayed_d;
X=convm(x,filter_length);
[A,E,y]= lms(X,d,miu,filter_length);

subplot(4,2,5)
plot(0:199,y(1:200),'g')
title('estimated signal y(t)');

subplot(4,2,7)
plot(0:199,E(1:200),'r')
title('error  e(t)');


subplot(4,2,6)
plot(abs( fftshift(fft(y) ) ))
title('filter output y(t)');


subplot(4,2,8)
plot(abs( fftshift(fft(E) ) ))
title('error');
