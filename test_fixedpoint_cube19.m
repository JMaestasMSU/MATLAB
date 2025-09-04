function [x, iter, p_estimate, errHist] = fixedPoint(g, x0, tol, maxiter)
% fixedPoint - Fixed-point iteration for solving x = g(x)
%
% This function implements the fixed-point iteration method to find roots.
% The method converts a root-finding problem f(x) = 0 into a fixed-point 
% problem x = g(x) by rearranging the original equation.
% Usage:
%   [x, iter, errHist] = fixedPoint(g, x0)
%   [x, iter, errHist] = fixedPoint(g, x0, tol)
%   [x, iter, errHist] = fixedPoint(g, x0, tol, maxiter)
% Inputs:
%   g           - function handle representing the iteration function g(x)
%                 Example: g = @(x) sqrt(10/(x+4)) for solving x^3+4x^2-10=0
%   x0          - starting value (scalar) - initial guess for the fixed point
%   tol         - tolerance for stopping (default 1e-6)
%                 Stopping criterion: |x_n - p_estimate| <= tol
%   maxiter     - maximum number of iterations allowed (default 100)
%              Safety limit to prevent infinite loops if method diverges
% Outputs:
%   x           - approximate fixed point (final iterate value)
%   iter        - number of iterations performed to reach convergence
%   errHist     - vector of estimated errors |x_n - p_estimate| at each iteration
%   p_estimate  - estimated fixed point p
% Method Overview:
%   1. Pre-estimate the fixed point by running a few iterations
%   2. Use this estimate to compute |p_n - p| style errors
%   3. Stop when estimated error falls below tolerance

    % Handle optional input arguments with default values (no tol sent)
    if nargin < 3 || isempty(tol)
        tol = 1e-6;        % Default tolerance
    end
    if nargin < 4 || isempty(maxiter)
        maxiter = 100;     % Default maximum iterations
    end
    
    % Since we don't know the true fixed point p, we run a few iterations
    % to get a reasonable estimate. This allows us to compute |p_n - p|
    
    x_est = x0;  % Start estimation from the same initial point
    
    % Run up to 10 iterations to estimate the fixed point location
    for i = 1:10
        x_new = g(x_est);  % Apply the fixed-point function
        
        % If consecutive estimates are very close, we have a good estimate
        if abs(x_new - x_est) < 1e-10
            fprintf('  Fixed point estimate converged after %d pre-iterations\n', i);
            break;
        end
        
        x_est = x_new;  % Update estimate for next iteration
    end
    
    p_estimate = x_est;  % Our best guess for the true fixed point p
    
    % Now perform the actual fixed-point iteration using our estimate
    % to compute textbook-style errors |p_n - p|
    
    % Pre-allocate error history array for efficiency
    errHist = zeros(1, maxiter);
    
    % Initialize iteration: start from original initial guess
    x_old = x0;
    
    % Main iteration loop: x_{n+1} = g(x_n)
    for iter = 1:maxiter
        x_new = g(x_old);  % Apply fixed-point function: x_{n+1} = g(x_n)
        
        % Calculate |p_n - p| using our estimated fixed point
        error_estimate = abs(x_new - p_estimate);
        errHist(iter) = error_estimate;  % Store for output
        
        % Detect divergence: if iteration produces NaN or Inf, method has failed
        if isnan(x_new) || isinf(x_new)
            errHist = errHist(1:iter);  % Trim unused array elements
            x = x_new;
            warning('fixedPoint:Divergence', 'Fixed-point iteration diverged at iteration %d', iter);
            return;
        end
        
        % Textbook stopping criterion: |p_n - p| ≤ tolerance
        % We use our estimated p since the true p is unknown
        if error_estimate <= tol
            errHist = errHist(1:iter);  % Trim unused array elements
            x = x_new;  % Return the converged value
            return;
        end
        
        x_old = x_new;  % Set up for next iteration: x_n becomes x_{n-1}
    end
    
    % If we reach this point, the method didn't converge within maxiter iterations
    errHist = errHist(1:maxiter);  % Keep full error history
    x = x_new;  % Return best available approximation
    warning('fixedPoint:NoConvergence', ...
            'Maximum iterations (%d) reached without meeting tolerance %.2e', ...
            maxiter, tol);
end