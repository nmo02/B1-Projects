function pi_estimate = MonteCarloPi(n_points)
    % Monte Carlo simulation to estimate the value of π
    inside_circle = 0;
    for i = 1:n_points
        x = rand;
        y = rand;
        if x^2 + y^2 <= 1
            inside_circle = inside_circle + 1;
        end
    end
    pi_estimate = (inside_circle / n_points) * 4;
end