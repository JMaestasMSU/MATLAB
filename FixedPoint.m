% Simple implementation of the fixed point root-finding method.

function [x, iter, errHist] = FixedPoint(g, x0, tol, maxiter)
% FixedPoint  Fixed-point iteration for solving x = g(x)
%
% Usage:
%   [x, iter, errHist] = FixedPoint(g, x0)
%   [x, iter, errHist] = FixedPoint(g, x0, tol)
%   [x, iter, errHist] = FixedPoint(g, x0, tol, maxiter)
%
% Inputs:
%   g        - function handle representing the iteration function g(x)
%   x0       - starting value (scalar)
%   tol      - tolerance for stopping (default 1e-6). Stopping criterion:
%              |x_{n+1} - x_n| <= tol
%   maxiter  - maximum number of iterations allowed (default 100)
% Outputs:
%   x        - approximate fixed point (last iterate)
%   iter     - number of iterations performed
%   errHist  - vector of absolute differences |x_{n+1}-x_n| at each iteration
%
% Notes / Comments:
%  - Fixed-point iteration computes x_{n+1} = g(x_n). Convergence is not
%    guaranteed for arbitrary g. A sufficient (but not necessary) condition
%    is |g'(x)| < 1 near the fixed point.
%  - This implementation uses the difference between successive iterates as
%    the stopping criterion because it is simple and commonly used.
%  - If you need a residual-based test, check |g(x)-x| after convergence.
%  - The function records errHist so you can inspect convergence behavior.

    % Provide default values for optional arguments
    if nargin < 3 || isempty(tol)
        tol = 1e-6;
    end
    if nargin < 4 || isempty(maxiter)
        maxiter = 100;
    end

    % Preallocate error history for efficiency (will trim if fewer iterations)
    errHist = zeros(1, maxiter);

    % Initialize
    x_old = x0;

    % Main fixed-point iteration loop:
    for iter = 1:maxiter
        % Compute next iterate
        x_new = g(x_old);

        % Compute difference between successive iterates (simple convergence measure)
        err = abs(x_new - x_old);
        errHist(iter) = err;

        % Check for invalid numbers (divergence / numerical blow-up)
        if isnan(err) || isinf(err)
            % Trim history and return with warning
            errHist = errHist(1:iter);
            x = x_new;
            warning('FixedPoint:NaNorInf', 'Iteration produced NaN or Inf at iter = %d', iter);
            return;
        end

        % Stopping criterion: successive iterates are within tolerance
        if err <= tol
            % Trim error history to actual length
            errHist = errHist(1:iter);
            x = x_new;
            return;
        end

        % Prepare for next iteration
        x_old = x_new;
    end

    % If we reach here, maximum iterations were used without meeting tol
    errHist = errHist(1:maxiter);
    x = x_new;
    warning('FixedPoint:NoConverge', 'Maximum iterations (%d) reached without meeting tolerance %g', maxiter, tol);
end