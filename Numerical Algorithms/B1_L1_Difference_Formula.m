
close all; clc


h=10.^[0:-0.5:-12];
exact_derivative=1;
forward_derivative=(exp(h)-1)./h;
backward_derivative=(1-exp(-h))./h;
central_derivative=(exp(h)-exp(-h))./(2*h);
difference_forward=abs(exact_derivative - forward_derivative);
difference_backward=abs(exact_derivative - backward_derivative);
difference_central=abs(exact_derivative - central_derivative);
figure
loglog(h,difference_forward,'-o',h,difference_backward,'-o',h, difference_central,'-o')
legend('forward','backward','central') 
xlabel('h')
ylabel('difference')
