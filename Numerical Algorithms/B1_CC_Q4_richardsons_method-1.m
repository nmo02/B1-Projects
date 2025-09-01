clear all;
close all;
clc;

% We want to estimate dy/dx at x=0.3679
start_point=0.3679;
end_point=1.0;

% To use Richardson extrapolation we need several values of step size
step=0.1;
max_step=0.5;
min_step=0.2;

% The actual function and derivative
x = 0:0.01:end_point;
y = (x.^x);
dydx = (x.^x).*(log(x)+1);
dydx_exact=0;

% Plot of f(x)
figure;
plot(x,y,'-');
legend('x^x'); 
xlabel('x'),ylabel('f(x)');
title('f(x)=x^x');

% Plot of f'(x)
figure;
plot(x,dydx,'-');
legend('x^x(log(x)+1)'); 
xlabel('x'),ylabel('df(x)/dx');
title('df(x)/dx=x^x(log(x)+1)');

% Arrays to hold different step size values and forward difference 
% values at those points.
forward_difference_size = [];
h_size = [];

% Open a figure to plot to
figure;

% Next step through the step sizes
for h=max_step:-step:min_step
    
    % Values of x and x+h
    x=start_point;
    x_plus_h=x+h;

    % Forward difference approximation
    forward_difference = ((x_plus_h^x_plus_h)-(x^x))/h;
    % Store this step size and value of forward difference
    forward_difference_size = [ forward_difference_size forward_difference ];
    h_size = [ h_size h];
    
    % Clear the figure
    clf;
    
    % Plot the actual function and the forward difference
    subplot(3,1,1);
    fplot(@(x) (x.^x), [0,1], 'b');
    xlabel('x'),ylabel('y');
    hold on;
    y=(x^x);
    y_h=(x_plus_h^x_plus_h);
    plot([x,x_plus_h],[y,y_h],'-o');
    legend('f(x)=x^x)'); 
    axis([0 1 0.6 1]);
    
    % Plot the actual df(x)/dx against the forward difference 
    % approximation at the x point of interest x=1/e
    subplot(3,1,2); 
    fplot(@(x) (x.^x).*(log(x)+1), [0,1], 'b');
    xlabel('x'),ylabel('y');
    hold on;
    plot(x,forward_difference,'-o');
    legend('df(x)/dx=x^x(log(x)+1))','Forward difference approximation'); 
    axis([0.365 0.375 -0.1 0.5]);

    % Plot the forward difference against step size
    subplot(3,1,3); 
    plot(h_size,forward_difference_size,'-o');
    axis([0.0 max_step+0.1 -0.05 0.5]);
    txt = ['Forward difference abs error = '...
        num2str(abs(dydx_exact - forward_difference))'.'];
    text(0.25,0.10,txt);

    pause(0.75);
    drawnow;
end

% Fit a polynomial (straight line for O(h) error in our algorithm 
P = polyfit(h_size,forward_difference_size,1);

% Plot the extrapolation based on the above fit
x=[0:0.001:h_size(1)];
yfit = P(1).*x+P(2);
P(2)
hold on;
plot(x,yfit,'r-.');

% Calculate the Richardson extrapolation error and print
r_error=abs(dydx_exact - P(2));
txt = ['Richardson extrapolation abs error = ' num2str(r_error)'.'];
text(0.025,0.415,txt);
drawnow;
