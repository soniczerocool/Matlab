     %System Identification of a 32-coefficient FIR filter 
     %(500 iterations).
     x  = randn(1,500);     % Input to the filter
     b  = fir1(31,0.5);     % FIR system to be identified
     n  = 0.1*randn(1,500); % Observation noise signal
     d  = filter(b,1,x)+n;  % Desired signal
     mu = 0.008;            % LMS step size
     h = adaptfilt.lms(32,mu);
     [y,e] = filter(h,x,d);
     subplot(2,1,1); plot(1:500,[d;y;e]);
     title('System Identification of an FIR filter');
     legend('Desired','Output','Error');
     xlabel('time index'); ylabel('signal value');
     subplot(2,1,2); stem([b.',h.Coefficients.']);
     legend('Actual','Estimated'); 
     xlabel('coefficient #'); ylabel('coefficient value'); grid on;