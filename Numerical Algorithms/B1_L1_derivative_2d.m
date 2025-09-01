close all; clc

lims=3;
figure;

for step=1:-0.05:0.05

    x=-lims:step:lims;
    y=-lims:step:lims;

    [X,Y]=meshgrid(x,y);

    gaussian=exp(-(X.*X+Y.*Y));

    XPlus=X+step;
    YPlus=Y+step;
    XMinus=X-step;
    YMinus=Y-step;

    difference=(exp(-(XPlus.*XPlus+YPlus.*YPlus))-exp(-(XPlus.*XPlus+YMinus.*YMinus))-exp(-(XMinus.*XMinus+YPlus.*YPlus))+exp(-(XMinus.*XMinus+YMinus.*YMinus)))/4*0.1*0.1;

    surfc(X,Y,difference)

    xlabel('X');
    ylabel('Y');
    zlabel('difference');
    
    drawnow
    pause(step);
    
end