clear all;
close all;
clc;

% Start and end points 
x_start=0.0;
x_end=11.0;

% The analytic solution to compare numerical results
xx=[x_start:0.001:x_end]';
fn=cos(xx);

% Let's work with a fixed number of steps, 40 is good.
numSteps=40;
    
% Calculate our step size
h=(x_end-x_start)/numSteps;
    
% The x points for us to evaluate at
x=[x_start:h:x_end]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Now we make two guesses dy(0)/dx = 1 and -1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First Guess...

% Size our vectors
y1_pc=x;
dydx1_pc=x;
   
% Initial Value y(0)=1
y1_pc(1)=1;

% Initial guess dy(0)/dx = 1
dydx1_pc=1;
  
% Now perform predictor corrector
for i=1:(length(x)-1)
    ypred=y1_pc(i)+h*dydx1_pc(i);
    dydxpred=dydx1_pc(i)-h*y1_pc(i);
    y1_pc(i+1)=y1_pc(i)+(h/2.0)*(dydx1_pc(i)+dydxpred);
    dydx1_pc(i+1)=dydx1_pc(i)-(h/2.0)*(y1_pc(i)+ypred);
end

% Plot our results for the first guess
figure;
plot(xx,fn);
hold on
plot(x,y1_pc,'o');
hold on
  
legend('Exact','y^{\prime}(0)=1');
xlabel('x');
ylabel('y');
axis([x_start x_end -3 3]);
txt = ['step size (h) = ' num2str(h)'.'];
text(0.2,-2.0,txt);
txt = ['number of steps = ' num2str(numSteps)'.'];
text(0.2,-2.3,txt);

% Second Guess...

% Size our vectors
y2_pc=x;
dydx2_pc=x;
   
% Initial Value y(0)=1
y2_pc(1)=1;

% Next guess dy(0)/dx = -1
dydx2_pc=-1;
  
% Now perform predictor corrector
for i=1:(length(x)-1)
    ypred=y2_pc(i)+h*dydx2_pc(i);
    dydxpred=dydx2_pc(i)-h*y2_pc(i);
    y2_pc(i+1)=y2_pc(i)+(h/2.0)*(dydx2_pc(i)+dydxpred);
    dydx2_pc(i+1)=dydx2_pc(i)-(h/2.0)*(y2_pc(i)+ypred);
end
 
% Plot our results for the second guess
figure;
plot(xx,fn);
hold on
plot(x,y2_pc,'o');
hold on
  
legend('Exact','y^{\prime}(0)=-1');
xlabel('x');
ylabel('y');
axis([x_start x_end -3 3]);
txt = ['step size (h) = ' num2str(h)'.'];
text(0.2,-2.0,txt);
txt = ['number of steps = ' num2str(numSteps)'.'];
text(0.2,-2.3,txt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Use our two guesses to interpolate a value for dy(0)/dx
% 
% ??? HOW ???
%
% 1. I know that y(11)=0 for the correct value of dy(0)/dx
%       -> It's one of my boundary conditions.
%
% 2. I know when dy(0)/dx =  1, y(11) = a
% 3. I know when dy(0)/dx = -1, y(11) = b
%
% 4. Now, "plot" y(11) values vs dy(0)/dx values, 
%       draw a straight line between them.
%
% 5. We know y(11)=0 for the correct value of dy(0)/dx
%       -> Read off the "correct" value of dy(0)/dx, where y(11) = 0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create the plot
figure;
ylabel('y^{\prime}(0) guess = {-1 or 1}, y=mx+c');
xlabel('y(11)=x');
hold on
   
% Find "correct" value of dy(0)/dx, where y(11) = 0
grad=(y1_pc(i+1) - y2_pc(i+1))/(dydx1_pc(1)-dydx2_pc(1));
const=y1_pc(i+1)-dydx1_pc(1)*grad;

% Plot the result
plot([0],[const],'o');
legend('y^{\prime}(0) guess = {-1 or 1}','Interpolated value for y^{\prime}(0)');
   
% Use the "correct/estimate" value of dy(0)/dx and perform PC again...

ye_pc=x;
dydxe_pc=x;
 
ye_pc(1)=1;
dydxe_pc=const;

for i=1:(length(x)-1)
    ypred=ye_pc(i)+h*dydxe_pc(i);
    ydpred=dydxe_pc(i)-h*ye_pc(i);
    ye_pc(i+1)=ye_pc(i)+(h/2.0)*(dydxe_pc(i)+ydpred);
    dydxe_pc(i+1)=dydxe_pc(i)-(h/2.0)*(ye_pc(i)+ypred);
end
   
% Finally plot all three solutions to compare
figure;
plot(xx,fn);
hold on
plot(x,y1_pc,'o');
hold on
plot(x,y2_pc,'+');
hold on
plot(x,ye_pc,'*');
  
legend('Exact','y^{\prime}(0)=1','y^{\prime}(0)=-1','y^{\prime}(0)=0.1434');
axis([x_start x_end -3 3]);
txt = ['step size (h) = ' num2str(h)'.'];
text(0.2,-2.0,txt);
txt = ['number of steps = ' num2str(numSteps)'.'];
text(0.2,-2.3,txt);
    
drawnow;