close all; clc

% Create a sine wave with some gaussian noise.
x = linspace(0, 4*pi, 20);
xVal = pi/10:(4*pi/19):4*pi+pi/10
y = sin(x) + 0.25*randn(size(x));
yVal = sin(xVal) + 0.25*randn(size(xVal));

% Create a vector to hold MSE
err = [];
errVal = [];

% A figure to plot to
hFig = figure(1);
set(hFig, 'Position', [400 100 500 1000])

% Itterate over polynomical fits 
for n=1:1:16
    
    % fit an nth degree polynomial to our noisey sine wave
    p=polyfit(x,y,n);

    % Evaluate the poly fit at the points of our data
    y1 = polyval(p,x);
    yTest = polyval(p,xVal);
    
    % Plot our noisey data
    subplot(3,1,1);
    plot(x,y,'*')
    hold on
    
    % Plot our model predictions
    plot(x,y1,'o')
    hold off
    ylim([-2,2]);
    title('Training data and model predictions')
    xlabel('x')
    ylabel('y')
    txt = sprintf('degree of polynomial\nn = %d', n);
    text(4.5,1.5,txt)
    legend('training', 'prediction');
    
    subplot(3,1,2);
    xm = linspace(0,4*pi);
    ym = polyval(p,xm);
    plot(xm,ym)
    hold on
    
    plot(x,y,'*')
    hold on
    
    hold on
    plot(xVal,yVal,'s');
    hold off
    ylim([-2,2]);
    title('Validation data and model predictions')
    xlabel('x')
    ylabel('y')
    txt = sprintf('degree of polynomial\nn = %d', n);
    text(4.5,1.5,txt)
    legend('prediction', 'training', 'validation');

    % Plot our MSE for each increase in degree of fit
    subplot(3,1,3); 
    err(end+1) = immse(y,y1);
    plot(err,'-+');
    hold on
     
    % Plot our MSE for each increase in degree of fit
    subplot(3,1,3); 
    errVal(end+1) = immse(yVal,yTest);
    plot(errVal,'-p');
    hold on
    %text(7,0.4,sprintf('degree of polynomial\nn = %d', n));
    xlim([0,16]);
    ylim([0, 0.6]);
    title('MSE vs nth degree of polynomial for training and validation data')
    xlabel('n')
    ylabel('MSE')
    legend('training', 'validation');

%     % Plot the real signal and the model
%     subplot(3,1,3);
%     yreal = sin(x1);
%     plot(x1,yreal)
%     hold on
%     xm = linspace(0,4*pi);
%     ym = polyval(p,xm);
%     plot(xm,ym)
%     hold off
%     ylim([-2,2]);
%     title('Real sine wave and model prediction')
%     xlabel('x')
%     ylabel('y')

    pause(0.5)
end

% figure
%  yreal = sin(x1);
%     plot(x1,yreal)
%     hold on
%     xm = linspace(0,4*pi);
%     ym = polyval(p,xm);
%     plot(xm,ym)
%     hold off
%     ylim([-2,2]);
%     title('Real sine wave and model prediction')
%     xlabel('x')
%     ylabel('y')
