function R=covar(x,p)
x=x(:);
m=lenght(x);
x=x-ones(m,1)*(sum(x)/m);
R=convm(x,p+1)'*convm(x,p+1)/(m-1);
end;