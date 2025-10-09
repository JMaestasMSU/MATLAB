% Define the integrand
f = @(x) x.^6 - x.^4;

% Integration bounds
a = 2;
b = 6;

% Calculate exact solution for comparison
% Antiderivative: x^7/7 - x^5/5
exact = (b^7/7 - b^5/5) - (a^7/7 - a^5/5);

fprintf('=== NUMERICAL INTEGRATION METHODS ===\n\n');
fprintf('Exact solution: %.10f\n\n', exact);

%% Method 1: Composite Trapezoidal Rule
fprintf('--- Composite Trapezoidal Rule ---\n');
n_trap = 100; % number of subintervals
h_trap = (b - a) / n_trap;

% Generate evaluation points
x_trap = linspace(a, b, n_trap + 1);
y_trap = f(x_trap);

% Composite trapezoidal: (h/2)[f(x0) + 2f(x1) + 2f(x2) + ... + 2f(xn-1) + f(xn)]
trap_result = (h_trap/2) * (y_trap(1) + 2*sum(y_trap(2:end-1)) + y_trap(end));

trap_error = abs(exact - trap_result);
fprintf('Result with n=%d: %.10f\n', n_trap, trap_result);
fprintf('Absolute error: %.10e\n', trap_error);
fprintf('Relative error: %.10e\n\n', trap_error/abs(exact));

%% Method 2: Composite Simpson's Rule
fprintf('--- Composite Simpson''s Rule ---\n');
n_simp = 100; % must be even
h_simp = (b - a) / n_simp;

% Generate evaluation points
x_simp = linspace(a, b, n_simp + 1);
y_simp = f(x_simp);

% Composite Simpson's: (h/3)[f(x0) + 4f(x1) + 2f(x2) + 4f(x3) + ... + f(xn)]
simp_result = (h_simp/3) * (y_simp(1) + 4*sum(y_simp(2:2:end-1)) + ...
                            2*sum(y_simp(3:2:end-2)) + y_simp(end));

simp_error = abs(exact - simp_result);
fprintf('Result with n=%d: %.10f\n', n_simp, simp_result);
fprintf('Absolute error: %.10e\n', simp_error);
fprintf('Relative error: %.10e\n\n', simp_error/abs(exact));

%% Method 3: Composite Midpoint Rule
fprintf('--- Composite Midpoint Rule ---\n');
n_mid = 100;
h_mid = (b - a) / n_mid;

% Evaluate at midpoints of each subinterval
midpoints = a + h_mid/2 + (0:n_mid-1)*h_mid;
y_mid = f(midpoints);

% Composite midpoint: h * sum(f(midpoint_i))
mid_result = h_mid * sum(y_mid);

mid_error = abs(exact - mid_result);
fprintf('Result with n=%d: %.10f\n', n_mid, mid_result);
fprintf('Absolute error: %.10e\n', mid_error);
fprintf('Relative error: %.10e\n\n', mid_error/abs(exact));

%% Method 4: Gaussian Quadrature (3-point)
fprintf('--- Gaussian Quadrature (3-point, n=3) ---\n');
fprintf('Using Gauss-Legendre quadrature on [-1,1]\n\n');

% Gauss-Legendre 3-point quadrature parameters
% Roots of Phi_3(x) = (1/2)(5x^3 - 3x)
fprintf('Step 1: Define nodes and weights on [-1,1]\n');
fprintf('  Nodes (roots of Phi_3): t_1 = -sqrt(3/5), t_2 = 0, t_3 = sqrt(3/5)\n');
fprintf('  Weights: a_1 = 5/9, a_2 = 8/9, a_3 = 5/9\n\n');

w = [5/9, 8/9, 5/9];
t_gauss = [-sqrt(3/5), 0, sqrt(3/5)];

fprintf('  t_1 = %.10f\n', t_gauss(1));
fprintf('  t_2 = %.10f\n', t_gauss(2));
fprintf('  t_3 = %.10f\n\n', t_gauss(3));

fprintf('  a_1 = %.10f\n', w(1));
fprintf('  a_2 = %.10f\n', w(2));
fprintf('  a_3 = %.10f\n\n', w(3));

% Transform integral from [a,b] to [-1,1]
% x = ((b-a)/2)*t + (a+b)/2, dx = ((b-a)/2)*dt
fprintf('Step 2: Map nodes from [-1,1] to [%d,%d]\n', a, b);
fprintf('  Transformation: x = ((b-a)/2)*t + (a+b)/2\n');
fprintf('  where (b-a)/2 = %.10f and (a+b)/2 = %.10f\n\n', (b-a)/2, (a+b)/2);

x_gauss = ((b-a)/2) * t_gauss + (a+b)/2;

fprintf('  x_1 = ((b-a)/2)*t_1 + (a+b)/2 = %.10f\n', x_gauss(1));
fprintf('  x_2 = ((b-a)/2)*t_2 + (a+b)/2 = %.10f\n', x_gauss(2));
fprintf('  x_3 = ((b-a)/2)*t_3 + (a+b)/2 = %.10f\n\n', x_gauss(3));

% Evaluate function at mapped nodes
fprintf('Step 3: Evaluate f(x) = x^6 - x^4 at each node\n');
y_gauss = f(x_gauss);

fprintf('  f(x_1) = f(%.10f) = %.10f\n', x_gauss(1), y_gauss(1));
fprintf('  f(x_2) = f(%.10f) = %.10f\n', x_gauss(2), y_gauss(2));
fprintf('  f(x_3) = f(%.10f) = %.10f\n\n', x_gauss(3), y_gauss(3));

% Compute weighted sum
fprintf('Step 4: Compute weighted sum\n');
fprintf('  Sum = a_1*f(x_1) + a_2*f(x_2) + a_3*f(x_3)\n');
weighted_sum = sum(w .* y_gauss);

fprintf('  Sum = (%.10f)*(%.10f) + (%.10f)*(%.10f) + (%.10f)*(%.10f)\n', ...
    w(1), y_gauss(1), w(2), y_gauss(2), w(3), y_gauss(3));
fprintf('  Sum = %.10f + %.10f + %.10f\n', ...
    w(1)*y_gauss(1), w(2)*y_gauss(2), w(3)*y_gauss(3));
fprintf('  Sum = %.10f\n\n', weighted_sum);

% Apply Jacobian factor
fprintf('Step 5: Apply Jacobian (b-a)/2 to account for transformation\n');
fprintf('  Integral = ((b-a)/2) * Sum\n');
gauss_result = ((b-a)/2) * weighted_sum;
fprintf('  Integral = (%.10f) * (%.10f)\n', (b-a)/2, weighted_sum);
fprintf('  Integral = %.10f\n\n', gauss_result);

gauss_error = abs(exact - gauss_result);
fprintf('Step 6: Compare with exact solution\n');
fprintf('  Exact result:    %.10f\n', exact);
fprintf('  Gaussian result: %.10f\n', gauss_result);
fprintf('  Absolute error:  %.10e\n', gauss_error);
fprintf('  Relative error:  %.10e\n\n', gauss_error/abs(exact));

fprintf('Note: 3-point Gaussian quadrature has degree of precision 2n-1 = 5\n');
fprintf('This means it integrates polynomials up to degree 5 exactly.\n');
fprintf('Our integrand (x^6 - x^4) is degree 6, so there is small error.\n\n');

%% Method 5: MATLAB's built-in integral function
fprintf('--- MATLAB Built-in integral() ---\n');
matlab_result = integral(f, a, b);
matlab_error = abs(exact - matlab_result);
fprintf('Result: %.10f\n', matlab_result);
fprintf('Absolute error: %.10e\n', matlab_error);
fprintf('Relative error: %.10e\n\n', matlab_error/abs(exact));

%% Comparison Summary
fprintf('=== SUMMARY OF METHODS ===\n');
fprintf('Method                      | Result         | Absolute Error\n');
fprintf('-----------------------------------------------------------\n');
fprintf('Exact                       | %14.10f | -\n', exact);
fprintf('Composite Trapezoidal (n=%d)| %14.10f | %.4e\n', n_trap, trap_result, trap_error);
fprintf('Composite Simpson (n=%d)    | %14.10f | %.4e\n', n_simp, simp_result, simp_error);
fprintf('Composite Midpoint (n=%d)   | %14.10f | %.4e\n', n_mid, mid_result, mid_error);
fprintf('Gaussian 3-point            | %14.10f | %.4e\n', gauss_result, gauss_error);
fprintf('MATLAB integral()           | %14.10f | %.4e\n', matlab_result, matlab_error);

%% Visualization
figure('Position', [100, 100, 1000, 600]);

% Plot the integrand
x_plot = linspace(a, b, 1000);
y_plot = f(x_plot);

subplot(2,1,1);
plot(x_plot, y_plot, 'b-', 'LineWidth', 2);
hold on;
area(x_plot, y_plot, 'FaceAlpha', 0.3);
xlabel('x');
ylabel('f(x) = x^6 - x^4');
title('Integrand: f(x) = x^6 - x^4');
grid on;
legend('f(x)', 'Area under curve');

% Plot error comparison
subplot(2,1,2);
methods = {'Trap.', 'Simp.', 'Mid.', 'Gauss', 'MATLAB'};
errors = [trap_error, simp_error, mid_error, gauss_error, matlab_error];
bar(errors);
set(gca, 'XTickLabel', methods);
ylabel('Absolute Error');
title('Error Comparison of Different Methods');
grid on;
set(gca, 'YScale', 'log');