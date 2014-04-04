%labtest2
%
clear all
f=10;
t=0:.001:10;
x=4*sin(2*pi*1*t);
d=awgn(x,10);
subplot(1,2,1)
plot(t,x)
subplot(1,2,2)
plot(t,d)

n0_delay=24;
h = dfilt.delay(n0_delay);
delayed_d=filter(h,d);
figure(2)
plot(t,delayed_d)



beta=1.5
mu=0.5
nord=30
p=nord
[A,E] = lms1(delayed_d,d,mu,nord)
figure(3)
t1=1:1:length(E);
plot(t1,E)
figure(4)
plot(t1,d)