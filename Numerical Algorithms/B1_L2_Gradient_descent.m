
close all; clc

learningRate = 0.01;

grad = 0.2;
intercept = 0;

N = 20;

% Create a sine wave with some gaussian noise.
x = linspace(0, 4*pi, N);
y = x.*grad + 0.25*randn(size(x)) + intercept;

m = 1.0;
b = mean(y);

mOld = m;
bOld = b;

finM = 1.0;
finB = 1.0;

cutoff = 0.00001;

n = 0;

figure;

while finM > cutoff && finB > cutoff
    mDer = 0;
    bDer = 0;
    for i=1:N
       
        % Calculate partial derivatives
        % -2x(y - (mx + b))
        mDer = mDer - (2.0*x(i)*(y(i) - (m*x(i) + b)));
    
        % -2(y - (mx + b))
        bDer = bDer - (2.0*(y(i) - (m*x(i) + b)));
    end

    % We subtract because the derivatives point in direction of steepest ascent
    m = m - (mDer / N) * learningRate;
    b = b - (bDer / N) * learningRate;
    
    finM = abs(m - mOld);
    finB = abs(b - bOld);
    
    mOld = m;
    bOld = b;
    
    ypred = m.*x + b;
    
    n = n + 1;
    
    plot(x,y,'*')
    hold on
    plot(x, ypred)
    ylim([-0.25 3.5]);
    hold off
    xlabel('x')
    ylabel('y')
    txt = sprintf('Itteration of Gradient Descent = %d', n);
    text(2.5,3.0,txt)
    legend('training', 'prediction');
    title('Linear regression using gradient descent')
    drawnow
    
end

txt = sprintf('Final m = %f and b = %f', m, b);
text(2.5,2.8,txt)
txt = sprintf('Errors in m = %f and b = %f', abs(m-grad), abs(b-intercept));
text(2.5,2.6,txt)