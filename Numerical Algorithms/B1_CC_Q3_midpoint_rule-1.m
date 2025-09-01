close all;
clear all;
clc;

% Integration range
start_point = 0.0;
end_point = 4.0;

% The exact analytic result
exact = 0.5;

% The gaussian function
xg = [start_point:.01:end_point];
%yg =normpdf(xg,0,1);
yg=(1.0/sqrt(2*pi))*exp(-((xg.*xg)/2.0));

% Number of steps to use
steps=10;

% Calculate our step size
h=(end_point-start_point)/steps;

% Create a vector that steps through our range
k=start_point:h:end_point;

% Set the sum that will hold our result to zero
int_mpr=0.0;

% Calculate the value of the function at each midpoint
y_point = exp(-((k+0.5*h).^2)/2)/sqrt(2*pi);

% Create a vector that holds the area of each strip
mpr=y_point.*h;

% Sum up all of the areas
int_mpr = sum(mpr, 'all')

% Calcluate the error between the exact analytic result and our midpoint
% rule approximation
error = (abs(int_mpr - exact)/exact).*100;
    
% Open a figure to plot to
figure;

% Create a fancy plot!
% Plot the rectangles
for i=1:1:steps
    lp=(i-1)*h+start_point;
    rp=(i)*h+start_point;
    pt=plot([lp rp], [y_point(i) y_point(i)], '--r');
    hold on
    plot([rp rp], [0 y_point(i)], '--r');
    hold on
    plot([lp lp], [0 y_point(i)], '--r');
    hold on
    txt = ['Step size (h) = ' num2str(h,'%4.4f') ];
    text(2.5,0.25,txt);
    txt = ['Percentage error = ' num2str(error,'%4.2f')];
    text(2.5,0.225,txt);
end

% Plot the actual function 
hold on
pe=plot(xg,yg,'-','LineWidth',2);
legend([pt pe],'Midpoint rule','Exact');
    
hold off
drawnow
