
mu=0.6;
gama=0.00001;
nord=12;
Fs=12000;
X=convm(x,nord);
[A,E,y]= llms(X,d,mu,gama,nord);


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
