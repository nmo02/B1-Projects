import numpy as np
import matplotlib.pyplot as plt

# Define Rosenbrock's function and its gradient and Hessian
def rosenbrock(x, y):
    return 100 * (y - x**2)**2 + (1 - x)**2

def rosenbrock_gradient(x, y):
    dfdx = -400 * x * (y - x**2) + 2 * (x - 1)
    dfdy = 200 * (y - x**2)
    return np.array([dfdx, dfdy])

def rosenbrock_hessian(x, y):
    d2fdx2 = 1200 * x**2 - 400 * y + 2
    d2fdy2 = 200
    d2fdxdy = -400 * x
    return np.array([[d2fdx2, d2fdxdy], [d2fdxdy, d2fdy2]])

# Define optimization methods
def gradient_descent(grad_func, start, learning_rate=0.001, max_iters=10000, tol=1e-6):
    x, y = start
    path = [(x, y)]
    for _ in range(max_iters):
        grad = grad_func(x, y)
        x -= learning_rate * grad[0]
        y -= learning_rate * grad[1]
        path.append((x, y))
        if np.linalg.norm(grad) < tol:
            break
    return np.array(path), len(path) - 1

def newtons_method(grad_func, hess_func, start, max_iters=6, tol=1e-6):
    x, y = start
    path = [(x, y)]
    for _ in range(max_iters):
        grad = grad_func(x, y)
        hess = hess_func(x, y)
        try:
            delta = np.linalg.solve(hess, grad)
        except np.linalg.LinAlgError:
            break
        x -= delta[0]
        y -= delta[1]
        path.append((x, y))
        if np.linalg.norm(grad) < tol:
            break
    return np.array(path), len(path) - 1

def gauss_newton(start, max_iters=5000, tol=1e-6):
    x, y = start
    path = [(x, y)]
    for _ in range(max_iters):
        r1 = 10 * (y - x**2)
        r2 = 1 - x
        r = np.array([r1, r2])
        
        dr1dx = -20 * x
        dr1dy = 10
        dr2dx = -1
        dr2dy = 0
        J = np.array([[dr1dx, dr1dy], [dr2dx, dr2dy]])
        
        try:
            delta = np.linalg.solve(J.T @ J, -J.T @ r)
        except np.linalg.LinAlgError:
            break
        
        x += delta[0]
        y += delta[1]
        path.append((x, y))
        
        if np.linalg.norm(delta) < tol:
            break
    return np.array(path), len(path) - 1

# Set initial points and grid for plotting
np.random.seed(42)  # For reproducible random points
initial_points = [np.random.uniform(-1, 1, 2) for _ in range(3)]
X, Y = np.meshgrid(np.linspace(-3, 3, 500), np.linspace(-3, 3, 500))
Z = rosenbrock(X, Y)

# Updated colors and plot settings
path_colors = ['blue', 'green', 'orange']
methods = [gradient_descent, newtons_method, gauss_newton]
method_names = ["Gradient Descent", "Newton's Method", "Gauss-Newton Method"]

# Create a single figure with three subplots
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle("Optimization of Rosenbrock's Function", fontsize=16)

for ax, method, name in zip(axes, methods, method_names):
    ax.set_facecolor("white")  # Set the background to white

    # Contour plot of the Rosenbrock function
    contour_filled = ax.contourf(X, Y, Z, levels=50, cmap="cividis", alpha=0.8)
    contour_lines = ax.contour(X, Y, Z, levels=15, colors="black", linewidths=0.5)
    ax.plot(1, 1, 'ro', markersize=10, label="Global Minimum (1,1)")

    # Plot paths for the current method and all starting points
    for color, start in zip(path_colors, initial_points):
        if name == "Newton's Method":
            path, iterations = method(rosenbrock_gradient, rosenbrock_hessian, start)
        elif name == "Gauss-Newton Method":
            path, iterations = method(start)
        else:
            path, iterations = method(rosenbrock_gradient, start)
        
        ax.plot(path[:, 0], path[:, 1], marker="o", markersize=3, color=color, label=f"Start {start}, Iterations: {iterations}")
        ax.arrow(path[0, 0], path[0, 1], path[1, 0] - path[0, 0], path[1, 1] - path[0, 1],
                 head_width=0.1, head_length=0.1, fc=color, ec=color)

    # Finalize individual subplot settings
    ax.set_xlim(-3, 3)
    ax.set_ylim(-3, 3)
    ax.set_title(f"{name}", fontsize=14)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend()

plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout to fit the title
plt.show()
