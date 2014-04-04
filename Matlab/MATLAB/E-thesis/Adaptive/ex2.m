
n0_delay=24;
h = dfilt.delay(n0_delay);
delayed_d=filter(h,d);

signal = sin(2*pi*24000*[0:12000-1]');
ha=adaptfilt.lms(12,0.006);
[y,e]=filter(ha,delayed_d,d);
% e(n) is the desired signal and it is not zero
plot(0:1999,y(1:2000),'b',0:1999,e(1:2000),'r',0:1999,signal(1:2000),'g')



