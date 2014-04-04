
% signal = sin(2*pi*24000*[0:12000-1]');
ha=adaptfilt.lms(12,0.006);
[y,e]=filter(ha,x,d);
% e(n) is the desired signal and it is not zero
plot(0:1999,y(1:2000),'b',0:1999,e(1:2000),'r')



%============================================

% 
% mu = .1;                % NLMS step size
% offset = .50;           % NLMS offset
% ha = adaptfilt.nlms(12,mu,offset);
% [y,e]=filter(ha,x,d);
% plot(0:1999,y(1:2000),'b',0:1999,e(1:2000),'r')
% 


%============================================















