% Simple implementation of the bisection root-finding method.
% mybisect(f,a,b,tol) returns c: an approximation to a root of f in [a,b]
% Preconditions: f(a)*f(b) < 0 (there is at least one sign change).

function c = mybisect(f, a, b, tol)
   % mybisect - find a root of f on [a,b] using the bisection method
   % Inputs:
   %   f   - function handle, f(x)
   %   a,b - interval endpoints (must satisfy f(a)*f(b) < 0)
   %   tol - tolerance for stopping criterion (half-interval width)
   % Output:
   %   c   - approximate root (midpoint of final interval)
   
   % Evaluate the function at the interval endpoints once.
   fa = f(a); % f at lower bound
   fb = f(b); % f at upper bound

   % Check the necessary condition for bisection: opposite signs at endpoints.
   % If not satisfied, the algorithm cannot guarantee a root in [a,b].
   if sign(fa)*sign(fb) >= 0
       error('f(a)f(b)<0 not satisfied!'); % signal misuse / bad interval
    end

   % Main bisection loop:
   % Continue until the half-width of the interval is <= tol.
   % The half-width (b-a)/2 is the maximum possible error of the midpoint.
   while (b - a)/2 > tol
       % Compute midpoint of the current interval
       c = (a+b)/2;
       fc = f(c); % evaluate f at midpoint

       % If fc is exactly zero (rare for floating point), we've found the root.
       if fc == 0
           return;
       end

       % Determine which subinterval contains the sign change:
       % If f(a) and f(c) have opposite signs, root is in [a,c], otherwise in [c,b].
       if sign(fa)*sign(fc) < 0
           % Root lies between a and c: update upper bound to c.
           b = c;
       else
           % Root lies between c and b: update lower bound to c.
           a = c;
           fa = fc; % update cached value for f(a)
       end
   end

   % When loop ends, return midpoint of final interval as the approximate root.
   c = (a+b)/2;
end

% Test function with roots inside the requested intervals:
% Tolerance requested by the assignment (accuracy within 10^-2)
tol = 1e-2;

% Define the function given in the problem:
f = @(x) x.^4 - 2.*x.^3 - 4.*x.^2 + 4.*x + 4;

% Intervals to test (each row is [a b])
intervals = [-2, -1; -1, 0; 0, 2; 2, 3];

% Loop over intervals, call mybisect, and display results
for k = 1:size(intervals,1)
    a = intervals(k,1);
    b = intervals(k,2);
    try
        % mybisect checks that f(a) and f(b) have opposite signs and then bisects
        root = mybisect(f, a, b, tol);
        fr = f(root);
        ok = abs(fr) <= tol; % check whether residual meets requested accuracy
        fprintf('Interval [%g, %g], root ≈ %.6f, f(root) = %.2e, |f|<=tol %d\n', a, b, root, fr, ok);
    catch E
        % If the interval does not bracket a root or another error occurs, print message
        fprintf('Interval [%g, %g]: ERROR - %s\n', a, b, E.message);
    end
end