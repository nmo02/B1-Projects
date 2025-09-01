% Define a list of sample sizes to study convergence
n_points_list = [10, 100, 1000, 10000, 100000, 1000000];
pi_estimates = zeros(size(n_points_list));
variances = zeros(size(n_points_list));

% Loop through each sample size to compute π estimates and variances
for i = 1:length(n_points_list)
    tic;
    [pi_estimates(i), variances(i)] = MonteCarloPiWithVariance(n_points_list(i));
    elapsed_time = toc;
    disp(['Monte Carlo time for ', num2str(n_points_list(i)), ' points: ', num2str(elapsed_time), ' seconds']);
end

% Calculate the error margins (95% confidence interval)
errors = 2 * sqrt(variances ./ n_points_list);  % 2*std_dev to approximate 95% CI

% Plot the convergence of π estimates with error bars
figure;
errorbar(n_points_list, pi_estimates, errors, '-o');
set(gca, 'XScale', 'log');  % Use a logarithmic scale for the x-axis
hold on;
yline(pi, '--r', 'True π'); % Horizontal line for the true value of π
xlabel('Number of Points');
ylabel('Estimate of π');
title('Convergence of Monte Carlo π Estimation with Error Bars');
legend('Monte Carlo Estimate', 'True π');
hold off;

% Display the final results
disp('Sample Sizes, π Estimates, and Variances:');
disp(table(n_points_list', pi_estimates', variances', errors', ...
    'VariableNames', {'SampleSize', 'PiEstimate', 'Variance', 'ErrorMargin'}));

% Local function for Monte Carlo π estimation with variance calculation
function [pi_estimate, variance] = MonteCarloPiWithVariance(n_points)
    inside_circle = 0;
    estimates = zeros(1, n_points);
    for i = 1:n_points
        x = rand;
        y = rand;
        if x^2 + y^2 <= 1
            inside_circle = inside_circle + 1;
        end
        estimates(i) = (inside_circle / i) * 4; % Update π estimate as we go
    end
    pi_estimate = estimates(end);  % Final π estimate
    variance = var(estimates);     % Calculate sample variance
end