% Clearing variables in memory and Matlab command screen
clear all;
close all;
clc;

% Dimensions of the simulation grid in x (xdim) and y (ydim) directions
xdim=20; ydim=100;

% Error cutoff
cutoff=0.001;

% Initializing previous (T_prev) and present (T_now) temp matrices
T_now  = zeros(xdim,ydim);
T_prev = zeros(xdim,ydim);

% Constant temperature of 67.5C applied to one boundary
T_now(1:xdim,1)  = 67.5;

% A simulation less ordinary 
%i=1:1:xdim; %x-co-ordinates for boundary 
%T_now(i,ydim)=exp((-1.*((xdim/2-i).*(xdim/2-i)))./(1.0.*xdim));
%
%j=1:1:ydim; %x-co-ordinates for boundary 
%T_now(xdim,j)=exp((-1.*((ydim/2-j).*(ydim/2-j)))./(2.0.*ydim));

% Ensure the current and previous grids are the same to start with. 
T_prev=T_now;

% Open a figure to plot to
figure;

% Set the iteration counter to one
t=1;

% Set a starting value for the error > cutoff;
error=2*cutoff;

% Keep iterating until the error between the current solution
% and last solution is below a predefined cutoff value
while error > cutoff
    
    % Solve for this time step, using the last time steps values 
    % over our solution space (mesh).
    for i=2:1:xdim-1
        for j=2:1:ydim-1
            T_now(i,j)=(T_prev(i+1,j)+T_prev(i-1,j)+T_prev(i,j+1)+T_prev(i,j-1))/4.0;
        end
    end
    
    % Calculate the difference between this time step and the last
    error=max(max(abs(T_now-T_prev)))
    
    % Now make this time step the last time step.
    T_prev=T_now;
    
    % Plot the current time step results
    surfc(T_now);
    title(sprintf('2D steady state heat equation at iteration = %1.2f',t),'Fontsize',11);
    xlabel('x','Fontsize',11); ylabel('y','Fontsize',11);
    zlabel(sprintf('T(x,y,iteration = %1.2f)',t),'Fontsize',11);
    drawnow;
    
    % Update the iteration counter and go again
    t=t+1;
end