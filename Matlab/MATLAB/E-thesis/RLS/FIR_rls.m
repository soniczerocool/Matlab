%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIR_rls : FIR RLS algorithm demo
% Author : Tamer Abdelazim Mellik
% Contact information : 
%Department of Electrical & Computer Engineering,
%University of Calgary,
%2500 University Drive N.W. ,
%Calgary, AB T2N 1N4 ,
%Canada .
% email :abdelasi@enel.ucalgary.ca  
% email : abdelazim@ieee.org
% Webpage : http://www.enel.ucalgary.ca/~abdelasi/
% Date    : 2-4-2002
% Updated : 30-10-2003
% Version : 1.1.0
% Reference : S. Haykin, Adaptive Filter Theory. 3rd edition, Upper Saddle River, NJ: Prentice-Hall, 1996. 
% Note : The author doesn't take any responsibility for any harm caused by the use of this file
clear all
close all
hold off
% Number of system points
N=2000;
inp = randn(N,1);
n = randn(N,1);
[b,a] = butter(2,0.25);
Gz = tf(b,a,-1);
%y(n)=b1*u(n)+b2*u(n-1)+b3*u(n-2)-a1*y(n-1)-a2*y(n-2)
h=[b -a(2:length(a))];


%If you don't have access to Control toolbox use the sample data
%load FIRsampledata;
%inp= FIRsampledata(1:2000);
%d= FIRsampledata(1:2000);

%channel system order fixed as we have 5 elements (3 in a and 2 in b)
inporder=3;
outorder=2;
sysorder = inporder + outorder ;

y = lsim(Gz,inp);
%add some noise
n = n * std(y)/(15*std(n));
d = y + n;
totallength=size(d,1);
%Take only 50 points for training ( N - inporder 47 = 50 - 3 )
N=50 ;	
%begin of the algorithm
%forgetting factor
lamda = 0.999 ;		
%initial P matrix
delta = 1e2 ;		 
P = delta * eye (sysorder ) ;
w = zeros ( sysorder  , 1 ) ;
for n = inporder : N 
    %u(n),u(n-1),u(n-2)
	u = inp(n:-1:n-inporder+1) ;
    %d(n-1),d(n-2)
    outp= d(n-1:-1:n-outorder) ;
    u=[u ; outp];
    phi = u' * P ;
	k = phi'/(lamda + phi * u );
    y(n)=w' * u;
    e(n) = d(n) - y(n) ;
	w = w + k * e(n) ;
	P = ( P - k * phi ) / lamda ;
    % Just for plotting
    Recordedw(1:sysorder,n)=w;
end 
%check of results
for n =  N+1 : totallength
    %u(n),u(n-1),u(n-2)
	u = inp(n:-1:n-inporder+1) ;
    %d(n-1),d(n-2)
    outp= d(n-1:-1:n-outorder) ;    
    u=[u ; outp];
    y(n) = w' * u ;
    e(n) = d(n) - y(n) ;
end 
hold on
plot(d)
plot(y,'r');
title('System output') ;
xlabel('Samples')
ylabel('True and estimated output')
figure
semilogy((abs(e))) ;
title('Error curve') ;
xlabel('Samples');
ylabel('Error value');
figure
plot(h, 'r+')
hold on
plot(w, '.')
legend('Filter weights','Estimated filter weights',4);
title('Comparison of the filter weights and estimated weights') ;
figure
plot(Recordedw(1:sysorder,sysorder:N)');
title('Estimated weights convergence') ;
xlabel('Samples');
ylabel('Weights value');
axis([1 N-sysorder min(min(Recordedw(1:sysorder,sysorder:N)')) max(max(Recordedw(1:sysorder,sysorder:N)')) ]);
hold off