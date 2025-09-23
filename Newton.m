function [c, iter, converged, history] = newton_improved(f, interval, varargin)
%NEWTON_IMPROVED Enhanced Newton-Raphson root finder with robust features
% [c, iter, converged, history] = NEWTON_IMPROVED(f, [a,b], options...)
%
% INPUTS:
%   f - function handle, scalar -> scalar
%   interval - [a,b], the interval containing the root
%   
% Optional name-value pairs:
%   'tol' - tolerance for |f(x)| [default 1e-8]
%   'step_tol' - tolerance for step size [default 1e-8]
%   'maxiter' - maximum iterations [default 100]
%   'df' - derivative function handle [default: numerical derivative]
%   'verbose' - display iteration progress [default false]
%   'initial_guess' - starting point [default: interval midpoint]
%
% OUTPUTS:
%   c - estimated root
%   iter - number of iterations used
%   converged - true if convergence achieved
%   history - struct with iteration history (x, f_vals, steps)

% Parse input arguments
p = inputParser;
addRequired(p, 'f', @(x) isa(x, 'function_handle'));
addRequired(p, 'interval', @(x) isnumeric(x) && length(x) == 2);
addParameter(p, 'tol', 1e-8, @(x) isnumeric(x) && x > 0);
addParameter(p, 'step_tol', 1e-8, @(x) isnumeric(x) && x > 0);
addParameter(p, 'maxiter', 100, @(x) isnumeric(x) && x > 0);
addParameter(p, 'df', [], @(x) isempty(x) || isa(x, 'function_handle'));
addParameter(p, 'verbose', false, @islogical);
addParameter(p, 'initial_guess', [], @(x) isempty(x) || isnumeric(x));

parse(p, f, interval, varargin{:});

tol = p.Results.tol;
step_tol = p.Results.step_tol;
maxiter = p.Results.maxiter;
df = p.Results.df;
verbose = p.Results.verbose;
initial_guess = p.Results.initial_guess;

% Validate interval
a = min(interval);
b = max(interval);

% Check Intermediate Value Theorem if possible
try
    fa = f(a);
    fb = f(b);
    if fa * fb > 0
        warning('IVT condition not satisfied: f(a)*f(b) > 0. Root may not exist in interval.');
    end
catch
    warning('Could not evaluate function at interval endpoints.');
end

% Set initial guess
if isempty(initial_guess)
    x = (a + b) / 2;  % midpoint
else
    x = initial_guess;
    if x < a || x > b
        warning('Initial guess outside interval. Using midpoint instead.');
        x = (a + b) / 2;
    end
end

% Initialize tracking variables
history.x = zeros(maxiter + 1, 1);
history.f_vals = zeros(maxiter + 1, 1);
history.steps = zeros(maxiter, 1);
converged = false;
iter = 0;

% Use analytical derivative if provided, otherwise numerical
if isempty(df)
    df = @(x) numerical_derivative(f, x);
end

if verbose
    fprintf('Newton-Raphson Method\n');
    fprintf('%-4s %-12s %-12s %-12s %-12s\n', 'Iter', 'x', 'f(x)', 'f''(x)', 'Step');
    fprintf('%s\n', repmat('-', 1, 60));
end

% Main Newton iteration loop
for iter = 1:maxiter
    try
        fval = f(x);
        dfval = df(x);
        
        % Store history
        history.x(iter) = x;
        history.f_vals(iter) = fval;
        
        if verbose
            fprintf('%-4d %-12.6f %-12.6e %-12.6e', iter-1, x, fval, dfval);
        end
        
        % Check convergence
        if abs(fval) < tol
            converged = true;
            if verbose
                fprintf(' %-12s\n', 'CONVERGED');
                fprintf('Converged due to |f(x)| < tol\n');
            end
            break;
        end
        
        % Check for problematic derivative
        if ~isfinite(dfval) || abs(dfval) < eps
            if verbose
                fprintf(' %-12s\n', 'BAD DERIV');
            end
            warning('Derivative is zero, infinite, or NaN at x = %.6f. Stopping.', x);
            break;
        end
        
        % Newton step
        step = -fval / dfval;
        x_new = x + step;
        
        history.steps(iter) = step;
        
        if verbose
            fprintf(' %-12.6e\n', step);
        end
        
        % Check for problematic new point
        if ~isfinite(x_new)
            warning('Next iterate is not finite. Stopping.');
            break;
        end
        
        % Check step size convergence
        if abs(step) < step_tol * (1 + abs(x))
            converged = true;
            x = x_new;
            if verbose
                fprintf('Converged due to small step size\n');
            end
            break;
        end
        
        x = x_new;
        
    catch ME
        warning('Error during iteration %d: %s', iter, ME.message);
        break;
    end
end

% Final function evaluation
if iter <= maxiter
    try
        fval = f(x);
        history.x(iter + 1) = x;
        history.f_vals(iter + 1) = fval;
    catch
        % If final evaluation fails, use last known value
    end
end

% Trim history arrays
history.x = history.x(1:iter+1);
history.f_vals = history.f_vals(1:iter+1);
history.steps = history.steps(1:iter);

c = x;

if ~converged && iter == maxiter
    warning('Maximum iterations (%d) reached without convergence.', maxiter);
end

if verbose
    fprintf('\nFinal result: x = %.10f, f(x) = %.6e, iterations = %d\n', c, f(c), iter);
end

end

function df_val = numerical_derivative(f, x)
%NUMERICAL_DERIVATIVE Compute numerical derivative using central difference
% Uses adaptive step size for better accuracy

h = max(1e-8, sqrt(eps) * (1 + abs(x)));

try
    df_val = (f(x + h) - f(x - h)) / (2 * h);
catch
    % Fallback to forward difference if central fails
    try
        df_val = (f(x + h) - f(x)) / h;
    catch
        df_val = NaN;
    end
end
end

%% Test Script
function run_newton_tests()
%RUN_NEWTON_TESTS Test the improved Newton's method on various problems

fprintf('Testing Improved Newton-Raphson Method\n');
fprintf('=====================================\n\n');

% Test problems
problems = {
    struct('name', 'i) x^3 - 2x^2 - 5 on [1,4]', ...
           'f', @(x) x.^3 - 2*x.^2 - 5, ...
           'df', @(x) 3*x.^2 - 4*x, ...
           'interval', [1, 4], ...
           'true_root', 2.690647448), % approximate
    
    struct('name', 'ii) x^3 + 3x^2 - 1 on [-3,-2]', ...
           'f', @(x) x.^3 + 3*x.^2 - 1, ...
           'df', @(x) 3*x.^2 + 6*x, ...
           'interval', [-3, -2], ...
           'true_root', -2.532088886), % approximate
    
    struct('name', 'iii) x - cos x on [0, π/2]', ...
           'f', @(x) x - cos(x), ...
           'df', @(x) 1 + sin(x), ...
           'interval', [0, pi/2], ...
           'true_root', 0.739085133), % approximate
    
    struct('name', 'iv) x - 0.8 - 0.2 sin x on [0, π/2]', ...
           'f', @(x) x - 0.8 - 0.2*sin(x), ...
           'df', @(x) 1 - 0.2*cos(x), ...
           'interval', [0, pi/2], ...
           'true_root', 1.066649016) % approximate
};

tol = 1e-10;  % High precision for testing
maxiter = 50;

fprintf('%-50s %-12s %-12s %-6s %-10s %-12s\n', ...
        'Problem', 'Root', 'Residual', 'Iters', 'Converged', 'Error');
fprintf('%s\n', repmat('-', 1, 100));

for i = 1:length(problems)
    prob = problems{i};
    
    % Test with analytical derivative
    [root_analytical, iter_analytical, conv_analytical] = ...
        newton_improved(prob.f, prob.interval, 'tol', tol, 'maxiter', maxiter, ...
                       'df', prob.df, 'verbose', false);
    
    % Test with numerical derivative
    [root_numerical, iter_numerical, conv_numerical] = ...
        newton_improved(prob.f, prob.interval, 'tol', tol, 'maxiter', maxiter, ...
                       'verbose', false);
    
    % Calculate errors if true root is known
    if isfield(prob, 'true_root')
        error_analytical = abs(root_analytical - prob.true_root);
        error_numerical = abs(root_numerical - prob.true_root);
    else
        error_analytical = NaN;
        error_numerical = NaN;
    end
    
    % Display results for analytical derivative
    fprintf('%-50s %-12.8f %-12.2e %-6d %-10s %-12.2e\n', ...
            [prob.name ' (analytical)'], root_analytical, abs(prob.f(root_analytical)), ...
            iter_analytical, char(string(conv_analytical)), error_analytical);
    
    % Display results for numerical derivative
    fprintf('%-50s %-12.8f %-12.2e %-6d %-10s %-12.2e\n', ...
            [prob.name ' (numerical)'], root_numerical, abs(prob.f(root_numerical)), ...
            iter_numerical, char(string(conv_numerical)), error_numerical);
    
    fprintf('\n');
end

% Demonstrate verbose output for one problem
fprintf('\nVerbose output example for problem i:\n');
fprintf('====================================\n');
prob = problems{1};
newton_improved(prob.f, prob.interval, 'tol', 1e-6, 'df', prob.df, 'verbose', true);

end

% Run tests if this file is executed directly
if ~exist('newton_test_suppress', 'var')
    run_newton_tests();
end