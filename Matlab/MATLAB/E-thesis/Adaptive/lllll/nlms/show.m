
beta=1;
nord=13;
X=convm(x,nord);
[A,E,y]= nlms(X,d,beta,nord);

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
