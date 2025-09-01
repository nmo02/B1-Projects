% Leibniz Series Method to estimate π and calculate error

% True value of π for comparison
true_pi = pi;

% Number of terms to use in the Leibniz series
n_terms = 100000;

% Start timing the Leibniz series calculation
tic;
% Calculate π using the Leibniz series
pi_estimate_leibniz = 0;
for k = 0:n_terms-1
    pi_estimate_leibniz = pi_estimate_leibniz + (-1)^k / (2 * k + 1);
end
pi_estimate_leibniz = pi_estimate_leibniz * 4;
leibniz_time = toc; % End timing

% Calculate the error by comparing with the true value of π
leibniz_error = abs(pi_estimate_leibniz - true_pi);

% Display the results
disp(['Leibniz series estimate of π with ', num2str(n_terms), ' terms: ', num2str(pi_estimate_leibniz)]);
disp(['Error for Leibniz series method: ', num2str(leibniz_error)]);
disp(['Leibniz series time: ', num2str(leibniz_time), ' seconds']);