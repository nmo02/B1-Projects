% Accurate root-finding script to estimate π using the bisection method

% True value of π for comparison
true_pi = pi;

% Define the function f(x) = sin(x), which has a root at π
f = @(x) sin(x);

% Set the interval and tolerance
a = 3;       % Lower bound of the interval
b = 4;       % Upper bound of the interval
tol = 1e-10; % High precision tolerance

% Check if f(a) and f(b) have opposite signs
if f(a) * f(b) > 0
    error('The function does not change sign over the interval. Choose a different interval.');
end

% Start timing the bisection method
tic;

% Bisection method to find the root near π
while (b - a) / 2 > tol
    midpoint = (a + b) / 2;
    f_mid = f(midpoint);
    
    if f_mid == 0
        root = midpoint;
        break;
    elseif f(a) * f_mid < 0
        b = midpoint; % Root is in the left subinterval
    else
        a = midpoint; % Root is in the right subinterval
    end
end
root = (a + b) / 2; % Final root estimate after loop

% Stop timing
root_finding_time = toc;

% Calculate the error for the root-finding method
root_finding_error = abs(root - true_pi);

% Display results
disp(['Root-finding estimate of π: ', num2str(root)]);
disp(['Error for root-finding method: ', num2str(root_finding_error)]);
disp(['Tolerance used: ', num2str(tol)]);
disp(['Root-finding time: ', num2str(root_finding_time), ' seconds']);