%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rlsalgo : IIR RLS algorithm demo
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

%If you don't have access to Control toolbox use the sample data
%load IIRsampledata;
%inp= IIRsampledata(1:2000);
%d= IIRsampledata(1:2000);

% use ldiv to get the approximate IIR weights of the filter ( a function only in the input)
% y=h*u
%ldiv is a function submitted to get inverse Z-transform (Matlab central file exchange)
%The first sysorder weight value
%use h=ldiv(b,a,sysorder)'; ==> here we use sysorder == 10

%channel system order (you can change the sysorder value and you don't need to change anything in the algorithm )
sysorder = 10 ;
%h= [0.0976   ; 0.2873  ;  0.3360   ; 0.2210   ; 0.0964   ; 0.0172 ;  -0.0159 ;  -0.0207  ; -0.0142  ; -0.0065 ; -0.0014 ;  0.0009 ;   0.0013   ; 0.0009   ; 0.0004  ;  0.0001  ; -0.0000  ; -0.0001  ; -0.0001  ; -0.0000];
h=[0.097631   0.287310   0.335965   0.220981 0.096354 0.017183  -0.015917 -0.020735  -0.014243  -0.006517 -0.001396   0.000856   0.001272  0.000914 0.000438 0.000108 -0.000044  -0.00008  -0.000058 -0.000029];
h=h(1:sysorder);
y = lsim(Gz,inp);
%add some noise
n = n * std(y)/(10*std(n));
d = y + n;
totallength=size(d,1);
%Take only 70 points for training ( N - systorder 70 = 80 - 10 )
N=80 ;	
%begin of the algorithm
%forgetting factor
lamda = 0.9995 ;		
%initial P matrix
delta = 1e10 ;		 
P = delta * eye (sysorder ) ;
w = zeros ( sysorder  , 1 ) ;
for n = sysorder : N 
	u = inp(n:-1:n-sysorder+1) ;
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
	u = inp(n:-1:n-sysorder+1) ;
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
legend('filter weights','Estimated filter weights');
title('Comparison of the filter weights and estimated weights') ;
figure
plot(Recordedw(1:sysorder,sysorder:N)');
title('Estimated weights convergence') ;
xlabel('Samples');
ylabel('Weights value');
axis([1 N-sysorder min(min(Recordedw(1:sysorder,sysorder:N)')) max(max(Recordedw(1:sysorder,sysorder:N)')) ]);
hold off