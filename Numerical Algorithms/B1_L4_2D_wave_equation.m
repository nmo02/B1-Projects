clc; clear all; close all;

% Phase velocity = 1
c = 1;  

% Step size
dx = 0.01;
dy = dx;

%Courant-Friedrich Stability Condition
sigma = 1/sqrt(2); 
gamma = 1/sqrt(2);

% time step
dt = sigma*(dx/c);

% number of points / mesh
t = 0:dt:1;
x = 0:dx:1;
y = 0:dy:1; 

% Set initial displacements to zero
u = zeros(length(x),length(y),length(t));

% Location of initial disturbance (stone in a pond for example)
xc=10.0;
yc=50.0;

% Varience for Gaussian disturbance
var=5;

% Initialisation of gaussian disturbance
for j=2:length(y)
    for i=2:length(x)
       exponent  = ((i-xc).^2 + (j-yc).^2)./(2*var^2);
        amplitude = 1 / (var * sqrt(2*pi));  
        val       = -0.2*amplitude  * exp(-exponent);
        u(i,j,2)  = val;
    end
end
 
% Begin time stepping
for n=2:length(t)-1
    for i=2:length(x)-1
        for j=2:length(y)-1
            u(i,j,n+1)= (sigma^2)*(u(i+1,j,n)-2*u(i,j,n)+u(i-1,j,n))...
                +(gamma^2)*(u(i,j+1,n)-2*u(i,j,n)+u(i,j-1,n)) + 2*u(i,j,n) - u(i,j,n-1);            
        end
    end
end

% Mesh for surfc plot
[X,Y] = meshgrid(x,y);

% Plot results
figure;
for j=1:length(t)
       p1 = surfc(X,Y,u(:,:,j)); 
       axis ([0 1 0 1 -0.05 0.035]);
       title(sprintf('2D wave equation at t = %1.2f',t(j)),'Fontsize',11);
       xlabel('x','Fontsize',11); ylabel('y','Fontsize',11);
       zlabel(sprintf('u(x,y,t = %1.2f)',t(j)),'Fontsize',11);
       %pause(1);
       drawnow;
end

