function df = numDiff(f, a, b, n)
% numDiff - Computes numerical derivatives using finite difference formulas
%
% Syntax: df = numDiff(f, a, b, n)
%
% Inputs:
%   f  - Function handle for the function to differentiate
%   a  - Left endpoint of interval [a,b]
%   b  - Right endpoint of interval [a,b]
%   n  - Number of subintervals (n+1 total points)
%
% Output:
%   df - Row vector of approximated derivative values at each grid point
%
% Method:
%   - Forward difference (O(h)) for first point: f'(x_0) = (f(x_0+h) - f(x_0))/h
%   - Centered difference (O(h^2)) for interior points: f'(x_i) = (f(x_i+h) - f(x_i-h))/(2h)
%   - Backward difference (O(h)) for last point: f'(x_n) = (f(x_n) - f(x_n-h))/h

    % Create equispaced grid
    h = (b - a) / n;              % Step size
    x = a:h:b;                     % Grid points from a to b (row vector)
    df = zeros(1, n+1);            % Initialize derivative vector (row vector)
    
    % Forward difference for first point (left endpoint)
    % Truncation error: O(h)
    df(1) = (feval(f, x(1) + h) - feval(f, x(1))) / h;
    
    % Centered difference for interior points
    % Truncation error: O(h^2) - significantly more accurate
    for i = 2:n
        df(i) = (feval(f, x(i) + h) - feval(f, x(i) - h)) / (2 * h);
    end
    
    % Backward difference for last point (right endpoint)
    % Truncation error: O(h)
    df(n+1) = (feval(f, x(n+1)) - feval(f, x(n+1) - h)) / h;
    
end

%% Define all functions and parameters
% Problem 1(a): f(x) = 2x^2 - 3x + 1 on [0,1]
f_a = @(x) 2*x.^2 - 3*x + 1;
df_a_exact = @(x) 4*x - 3;
a_a = 0; b_a = 1;

% Problem 1(b): f(x) = e^x - x on [0,1]
f_b = @(x) exp(x) - x;
df_b_exact = @(x) exp(x) - 1;
a_b = 0; b_b = 1;

% Problem 1(c): f(x) = e^(2x) on [1,2]
f_c = @(x) exp(2*x);
df_c_exact = @(x) 2*exp(2*x);
a_c = 1; b_c = 2;

% Problem 1(d): f(x) = e^(-4x) on [-1,0]
f_d = @(x) exp(-4*x);
df_d_exact = @(x) -4*exp(-4*x);
a_d = -1; b_d = 0;

% Test parameters
h_values = [0.1, 0.05, 0.01];

%% Storage structures for results
results = struct();

%% Problem 1(a): f(x) = 2x^2 - 3x + 1 on [0,1] - Compute and store
% Initialize storage arrays for this problem
results.a.h_vals = h_values;                      % Store the h values tested
results.a.max_errors = zeros(size(h_values));     % Maximum error at each h
results.a.ratios = zeros(size(h_values));         % Error ratios (to verify O(h²))
results.a.expected_ratios = zeros(size(h_values)); % Theoretical ratios

for i = 1:length(h_values)
    h = h_values(i);                    % Current step size
    n = round((b_a - a_a) / h);         % Number of subintervals
    x = a_a:h:b_a;                      % Grid points (row vector)
    
    % Compute numerical derivative using our numDiff function
    df_numerical = numDiff(f_a, a_a, b_a, n);
    
    % Compute exact derivative at grid points (ensure row vector)
    df_exact = df_a_exact(x);
    df_exact = df_exact(:)';  % Force to row vector
    
    % Calculate pointwise absolute errors
    errors = abs(df_numerical - df_exact);
    results.a.max_errors(i) = max(errors);
    
    % Compute error ratio to verify O(h²) convergence
    % When h is halved, error should decrease by factor of 4
    if i > 1
        results.a.ratios(i) = results.a.max_errors(i-1) / results.a.max_errors(i);
        results.a.expected_ratios(i) = (h_values(i-1) / h_values(i))^2;
    end
end

%% Problem 1(b): f(x) = e^x - x on [0,1] - Compute and store
% Initialize storage arrays
results.b.h_vals = h_values;
results.b.max_errors = zeros(size(h_values));
results.b.ratios = zeros(size(h_values));
results.b.expected_ratios = zeros(size(h_values));

for i = 1:length(h_values)
    h = h_values(i);
    n = round((b_b - a_b) / h);
    x = a_b:h:b_b;
    
    % Compute derivatives
    df_numerical = numDiff(f_b, a_b, b_b, n);
    df_exact = df_b_exact(x);
    df_exact = df_exact(:)';  % Force to row vector
    
    % Store maximum error
    errors = abs(df_numerical - df_exact);
    results.b.max_errors(i) = max(errors);
    
    % Compute convergence ratios
    if i > 1
        results.b.ratios(i) = results.b.max_errors(i-1) / results.b.max_errors(i);
        results.b.expected_ratios(i) = (h_values(i-1) / h_values(i))^2;
    end
end

%% Problem 1(c): f(x) = e^(2x) on [1,2] - Compute and store
% Initialize storage arrays
results.c.h_vals = h_values;
results.c.max_errors = zeros(size(h_values));
results.c.ratios = zeros(size(h_values));
results.c.expected_ratios = zeros(size(h_values));

for i = 1:length(h_values)
    h = h_values(i);
    n = round((b_c - a_c) / h);
    x = a_c:h:b_c;
    
    % Compute derivatives
    df_numerical = numDiff(f_c, a_c, b_c, n);
    df_exact = df_c_exact(x);
    df_exact = df_exact(:)';  % Force to row vector
    
    % Store maximum error
    errors = abs(df_numerical - df_exact);
    results.c.max_errors(i) = max(errors);
    
    % Compute convergence ratios
    if i > 1
        results.c.ratios(i) = results.c.max_errors(i-1) / results.c.max_errors(i);
        results.c.expected_ratios(i) = (h_values(i-1) / h_values(i))^2;
    end
end

%% Problem 1(d): f(x) = e^(-4x) on [-1,0] - Compute and store
% Initialize storage arrays
results.d.h_vals = h_values;
results.d.max_errors = zeros(size(h_values));
results.d.ratios = zeros(size(h_values));
results.d.expected_ratios = zeros(size(h_values));

for i = 1:length(h_values)
    h = h_values(i);
    n = round((b_d - a_d) / h);
    x = a_d:h:b_d;
    
    % Compute derivatives
    df_numerical = numDiff(f_d, a_d, b_d, n);
    df_exact = df_d_exact(x);
    df_exact = df_exact(:)';  % Force to row vector
    
    % Store maximum error
    errors = abs(df_numerical - df_exact);
    results.d.max_errors(i) = max(errors);
    
    % Compute convergence ratios
    if i > 1
        results.d.ratios(i) = results.d.max_errors(i-1) / results.d.max_errors(i);
        results.d.expected_ratios(i) = (h_values(i-1) / h_values(i))^2;
    end
end

%% Large h exploration (h >= 1) - Compute and store
% Test what happens when h is too large (violates h << 1 assumption)
% Use extended domain [0,4] so we have enough room for large steps
a_large = 0; b_large = 4;
large_h_values = [1.0, 1.5, 2.0];

% Initialize storage
results.large.h_vals = large_h_values;
results.large.max_errors = zeros(size(large_h_values));
results.large.n_vals = zeros(size(large_h_values));     % Track number of points
results.large.valid = true(size(large_h_values));       % Flag if enough points exist

for i = 1:length(large_h_values)
    h = large_h_values(i);
    n = round((b_large - a_large) / h);  % Number of subintervals
    results.large.n_vals(i) = n;
    
    % Check if we have enough points to apply centered differences
    % Need at least 2 subintervals (3 points) for the method to work
    if n < 2
        results.large.valid(i) = false;
        continue;  % Skip this case
    end
    
    % Compute derivatives on extended domain
    df_numerical = numDiff(f_b, a_large, b_large, n);  % Use f_b = e^x - x
    h = (b_large - a_large) / n;   % Use same step size as numDiff
    x = a_large:h:b_large;         % Create matching grid
    fprintf('Size of df_numerical: [%s]\n', num2str(size(df_numerical)));
    df_exact = df_b_exact(x);
    fprintf('Size of df_exact before reshape: [%s]\n', num2str(size(df_exact)));
    df_exact = df_exact(:)';
    fprintf('Size of df_exact after reshape: [%s]\n', num2str(size(df_exact)));
    
    % Store maximum error
    errors = abs(df_numerical - df_exact);
    results.large.max_errors(i) = max(errors);
end

%% Visualization
figure('Position', [100, 100, 1200, 800]);

% Plot 1(a)
subplot(2,2,1);
x = linspace(0, 1, 100);
plot(x, df_a_exact(x), 'b-', 'LineWidth', 2); hold on;
h = 0.1; n = (1-0)/h; x_num = 0:h:1;
plot(x_num, numDiff(f_a, 0, 1, n), 'ro-', 'MarkerSize', 6);
xlabel('x'); ylabel("f'(x)"); title('(a) f(x) = 2x² - 3x + 1');
legend('Exact', 'Numerical (h=0.1)');
grid on;

% Plot 1(b)
subplot(2,2,2);
x = linspace(0, 1, 100);
plot(x, df_b_exact(x), 'b-', 'LineWidth', 2); hold on;
h = 0.1; n = (1-0)/h; x_num = 0:h:1;
plot(x_num, numDiff(f_b, 0, 1, n), 'ro-', 'MarkerSize', 6);
xlabel('x'); ylabel("f'(x)"); title('(b) f(x) = e^x - x');
legend('Exact', 'Numerical (h=0.1)');
grid on;

% Plot 1(c)
subplot(2,2,3);
x = linspace(1, 2, 100);
plot(x, df_c_exact(x), 'b-', 'LineWidth', 2); hold on;
h = 0.1; n = (2-1)/h; x_num = 1:h:2;
plot(x_num, numDiff(f_c, 1, 2, n), 'ro-', 'MarkerSize', 6);
xlabel('x'); ylabel("f'(x)"); title('(c) f(x) = e^{2x}');
legend('Exact', 'Numerical (h=0.1)');
grid on;

% Plot 1(d)
subplot(2,2,4);
x = linspace(-1, 0, 100);
plot(x, df_d_exact(x), 'b-', 'LineWidth', 2); hold on;
h = 0.1; n = (0-(-1))/h; x_num = -1:h:0;
plot(x_num, numDiff(f_d, -1, 0, n), 'ro-', 'MarkerSize', 6);
xlabel('x'); ylabel("f'(x)"); title('(d) f(x) = e^{-4x}');
legend('Exact', 'Numerical (h=0.1)');
grid on;

%% Now print all results
fprintf('\n');
fprintf('========================================================\n');
fprintf('    NUMERICAL DIFFERENTIATION ANALYSIS RESULTS\n');
fprintf('========================================================\n\n');

%% Print Problem 1(a) results
fprintf('=== Problem 1(a): f(x) = 2x² - 3x + 1 on [0,1] ===\n\n');
for i = 1:length(results.a.h_vals)
    fprintf('h = %.4f:\n', results.a.h_vals(i));
    fprintf('  Maximum error: %.8e\n', results.a.max_errors(i));
    if i > 1
        fprintf('  Error ratio: %.4f (expected ~%.4f for O(h²))\n', ...
                results.a.ratios(i), results.a.expected_ratios(i));
    else
        fprintf('  Error ratio: N/A (first case)\n');
    end
    fprintf('\n');
end

%% Print Problem 1(b) results
fprintf('=== Problem 1(b): f(x) = e^x - x on [0,1] ===\n\n');
for i = 1:length(results.b.h_vals)
    fprintf('h = %.4f:\n', results.b.h_vals(i));
    fprintf('  Maximum error: %.8e\n', results.b.max_errors(i));
    if i > 1
        fprintf('  Error ratio: %.4f (expected ~%.4f for O(h²))\n', ...
                results.b.ratios(i), results.b.expected_ratios(i));
    else
        fprintf('  Error ratio: N/A (first case)\n');
    end
    fprintf('\n');
end

%% Print Problem 1(c) results
fprintf('=== Problem 1(c): f(x) = e^(2x) on [1,2] ===\n\n');
for i = 1:length(results.c.h_vals)
    fprintf('h = %.4f:\n', results.c.h_vals(i));
    fprintf('  Maximum error: %.8e\n', results.c.max_errors(i));
    if i > 1
        fprintf('  Error ratio: %.4f (expected ~%.4f for O(h²))\n', ...
                results.c.ratios(i), results.c.expected_ratios(i));
    else
        fprintf('  Error ratio: N/A (first case)\n');
    end
    fprintf('\n');
end

%% Print Problem 1(d) results
fprintf('=== Problem 1(d): f(x) = e^(-4x) on [-1,0] ===\n\n');
for i = 1:length(results.d.h_vals)
    fprintf('h = %.4f:\n', results.d.h_vals(i));
    fprintf('  Maximum error: %.8e\n', results.d.max_errors(i));
    if i > 1
        fprintf('  Error ratio: %.4f (expected ~%.4f for O(h²))\n', ...
                results.d.ratios(i), results.d.expected_ratios(i));
    else
        fprintf('  Error ratio: N/A (first case)\n');
    end
    fprintf('\n');
end

%% Print large h exploration results
fprintf('=== Exploration of large h values (h ≥ 1) on [0,4] ===\n\n');
fprintf('Testing f(x) = e^x - x:\n');
for i = 1:length(results.large.h_vals)
    h = results.large.h_vals(i);
    n = results.large.n_vals(i);
    
    if ~results.large.valid(i)
        fprintf('h = %.2f: Too few points (n = %d)\n', h, n);
        continue;
    end
    
    fprintf('h = %.2f (n = %d):\n', h, n);
    fprintf('  Maximum error: %.8e\n', results.large.max_errors(i));
    fprintf('  This is PROBLEMATIC because:\n');
    fprintf('    - Taylor series assumes h << 1\n');
    fprintf('    - Truncation error O(h²) becomes very large\n');
    fprintf('    - Lost local approximation property\n\n');
end

fprintf('========================================================\n');
fprintf('Analysis complete! Check the figure for visualizations.\n');
fprintf('========================================================\n');
legend('Exact', 'Numerical (h=0.1)');
grid on;