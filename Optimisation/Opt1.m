% Task 1: Visualization of Rosenbrock's Function

% Define the function
rosenbrock = @(x, y) 100 * (y - x.^2).^2 + (1 - x).^2;

% Generate mesh grid for x and y in the interval [-2, 2]
x = linspace(-2, 2, 400);
y = linspace(-2, 2, 400);
[X, Y] = meshgrid(x, y);

% Compute function values
Z = rosenbrock(X, Y);

% 3D Plot of Rosenbrock's Function
figure;
subplot(1, 2, 1);
surf(X, Y, Z, 'EdgeColor', 'none');
colormap('parula');
xlabel('X');
ylabel('Y');
zlabel('f(X, Y)');
title('3D Plot of Rosenbrock''s Function');
view(135, 30); % Adjust the viewing angle

% 2D Contour Plot of Rosenbrock's Function
subplot(1, 2, 2);
contourf(X, Y, Z, 50, 'LineColor', 'none');
colormap('parula');
colorbar;
xlabel('X');
ylabel('Y');
title('2D Contour Plot of Rosenbrock''s Function');
hold on;
% Mark the global minimum at (1, 1)
plot(1, 1, 'ro', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Minimum (1, 1)');
legend;
hold off;