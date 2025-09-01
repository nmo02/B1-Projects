clear all;
close all;
clc;

% Start and end points.
x_start=0.0;
x_end=11.0;

% The analytic solution to compare numerical results
xx=[x_start:0.001:x_end]';
fn=cos(xx);

% Open a figure to plot to
figure;

% Loop over a different number of steps to see how
% step size influences solution
for numSteps=2:1:100;
    
    % clear the figure
    clf;
    
    % Calculate step size
    h=(x_end-x_start)/numSteps;
    
    % The x range for our step size    
    x=[x_start:h:x_end]';
    
    % Size vectors for Euler and PC methods.
    y=x;
    dydx=x;
    
    y_pc=x;
    dydx_pc=x;

    % Initial values (as in lecture).
    y(1)=1;
    dydx=0;
    
    y_pc(1)=1;
    dydx_pc=0;
 
    % Perform Euler
    for i=1:(length(x)-1)
        y(i+1)=y(i)+h*dydx(i);
        dydx(i+1)=dydx(i)-h*y(i);
    end
    
    % Perform PC
    for i=1:(length(x)-1)
        ypred=y_pc(i)+h*dydx_pc(i);
        dydxpred=dydx_pc(i)-h*y_pc(i);
        
        y_pc(i+1)=y_pc(i)+(h/2.0)*(dydx_pc(i)+dydxpred);
        dydx_pc(i+1)=dydx_pc(i)-(h/2.0)*(y_pc(i)+ypred);
    end
 
    % Plot Euler and PC results for this step 
    plot(xx,fn);
    hold on
    plot(x,y,'*');        
    hold on
    plot(x,y_pc,'o');
    hold on
  
    % Axis, legend, step size adn number of steps 
    legend('Exact','Euler','Predictor-Corrector','Location','northwest');
    axis([x_start x_end -3 3]);
    txt = ['step size (h) = ' num2str(h)'.'];
    text(0.2,-2.0,txt);
    txt = ['number of steps = ' num2str(numSteps)'.'];
    text(0.2,-2.3,txt);
    
    drawnow;
    
    if mod(numSteps,20)==0
    %    pause(5);
    end
end