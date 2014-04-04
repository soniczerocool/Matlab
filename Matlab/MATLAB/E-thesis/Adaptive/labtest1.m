beta=1.5
mu=0.5
nord=30
p=nord
[A,E] = lms1(x,d,mu,nord)
figure(1)
t=1:1:length(E);
plot(t,E)
figure(2)
plot(t,d)