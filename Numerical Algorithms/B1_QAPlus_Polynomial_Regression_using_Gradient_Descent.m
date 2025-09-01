clear all;
close all;
clc

% Learning rate
learningRate = 0.00125    ;

% Polynomial degree
M = 3;

% Add one on to M for the intercept coefficient
M = M + 1;

% Polynomial coefficients
%alpha = normrnd(1,3,[1,M]);
alpha = [-2 -1 2 1];

% Number of training data
N = 40;

% Create a polynomial wave with some gaussian noise.
x = linspace(-4, 4, N);
y = polyval(alpha, x);

% Add some noise
y = y + (0.5)*M*randn(size(x));

delta_mse = 1;
last_mse = 1;
cutoff = 0.001;
n = 0;

alpha_pred = 1.5*randn(1,M);
alpha_pred = [0 0 0 1];

figure;

while delta_mse > cutoff
    alpha_der = zeros(1,M);
    for i=1:N 
       for j=1:M
            % Calculate partial derivatives
            % -2x^j(y - y_pred)
            alpha_der(j) = alpha_der(j) - (2.0*(x(i)^(M-j))*(y(i) - polyval(alpha_pred, x(i))));
        end
    end

    % We subtract because the derivatives point in direction of steepest ascent
    alpha_pred = alpha_pred - (alpha_der / N) * learningRate;
    
    ypred = polyval(alpha_pred, x);
      
    mse=immse(y,ypred)
    delta_mse = abs(mse - last_mse);
    last_mse = mse; 
    n = n + 1;
    
    plot(x,y,'*')
    hold on
    plot(x, ypred)
    hold off
    xlabel('x')
    ylabel('y')
    txt = sprintf('Itteration of Gradient Descent = %d', n);
    text(-3,125,txt)
    txt = sprintf('MSE error = %f', last_mse);
    text(-3,100,txt)
    legend('training', 'prediction');
    title('Polynomial regression using gradient descent')
    drawnow
    
end
