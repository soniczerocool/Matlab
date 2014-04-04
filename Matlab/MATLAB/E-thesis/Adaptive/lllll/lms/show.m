
mu=0.0361;
nord=22;
Fs=12000;

X=convm(x,nord);
[A,E,y]= lms(X,d,mu,nord);


subplot(4,1,1)
plot(x)
Title('x(n)')

subplot(4,1,2)
plot(y)
Title('y(n)')


subplot(4,1,3)
plot(d)
Title('d(n)')
subplot(4,1,4)
plot(E,'r')
Title('Desired Signal')
