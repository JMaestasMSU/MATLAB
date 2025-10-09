
function g = cubic_spline(x, f, bc, xt)
% cubic_spline - Compute cubic spline interpolation
%
% Inputs:
% x  - vector of x-coordinates (knot points)
% f  - vector of function values at x
% bc - boundary condition vector [val1, val2, type] where type: NaN=natural, 1=clamped
% xt - vector of evaluation points
%
% Outputs:
% g  - interpolated values at xt points

    % Number of spline segments
    n = length(x) - 1;
    
    % Initialize coefficient matrix A and right-hand side vector B
    % Each segment has 4 coefficients, so we have 4n unknowns
    A = zeros(4*n, 4*n);
    B = zeros(4*n, 1);

    % Build the linear system for each spline segment
    for i = 1:n
        % Calculate interval width
        h = x(i+1) - x(i);
        
        if i == n  % Last segment - handle boundary conditions
            % Constraint 1: S_i(x_i) = f_i (function value at left endpoint)
            A(4*i-3, 4*i-3:4*i) = [1, 0, 0, 0];
            B(4*i-3) = f(i);
            
            % Constraint 2: S_i(x_{i+1}) = f_{i+1} (function value at right endpoint)
            A(4*i-2, 4*i-3:4*i) = [1, h, h^2, h^3];
            B(4*i-2) = f(i+1);
            
            % Apply boundary conditions
            if isnan(bc(3))  % Natural boundary conditions (second derivative = 0)
                % Left boundary: S''_1(x_1) = 0
                A(4*i-1, 1:4) = [0, 0, 2, 0];
                B(4*i-1) = 0;
                
                % Right boundary: S''_n(x_{n+1}) = 0
                A(4*i, 4*i-3:4*i) = [0, 0, 2, 6*h];
                B(4*i) = 0;
                
            elseif bc(3) == 1  % Clamped boundary conditions (specified first derivatives)
                % Left boundary: S'_1(x_1) = bc(1)
                A(4*i-1, 1:4) = [0, 1, 0, 0];
                B(4*i-1) = bc(1);
                
                % Right boundary: S'_n(x_{n+1}) = bc(2)
                A(4*i, 4*i-3:4*i) = [0, 1, 2*h, 3*h^2];
                B(4*i) = bc(2);
                
            else
                % Other boundary condition types can be added here
            end
            
        else  % Interior segments - apply continuity constraints
            % Constraint 1: S_i(x_i) = f_i (function value at left endpoint)
            A(4*i-3, 4*i-3:4*i) = [1, 0, 0, 0];
            B(4*i-3) = f(i);
            
            % Constraint 2: S_i(x_{i+1}) = f_{i+1} (function value at right endpoint)
            A(4*i-2, 4*i-3:4*i) = [1, h, h^2, h^3];
            B(4*i-2) = f(i+1);
            
            % Constraint 3: First derivative continuity S'_i(x_{i+1}) = S'_{i+1}(x_{i+1})
            A(4*i-1, 4*i-3:4*i+4) = [0, 1, 2*h, 3*h^2, 0, -1, 0, 0];
            B(4*i-1) = 0;
            
            % Constraint 4: Second derivative continuity S''_i(x_{i+1}) = S''_{i+1}(x_{i+1})
            A(4*i, 4*i-3:4*i+4) = [0, 0, 2, 6*h, 0, 0, -2, 0];
            B(4*i) = 0;
        end
    end
    
    % Solve the linear system for spline coefficients
    C = A \ B;
    
    % Evaluate the spline at the query points
    m = length(xt);
    g = zeros(1, m);  % Pre-allocate output array
    
    for i = 1:m
        % Find which segment contains the evaluation point
        for j = 1:n
            if j == n  % Last segment - include right endpoint
                if xt(i) >= x(j) && xt(i) <= x(j+1)
                    % Calculate offset from left endpoint of segment
                    xh = xt(i) - x(j);
                    
                    % Evaluate cubic polynomial: a + b*xh + c*xh^2 + d*xh^3
                    g(i) = C(4*j-3) + C(4*j-2)*xh + C(4*j-1)*xh^2 + C(4*j)*xh^3;
                    break;
                end
            else  % Interior segments - exclude right endpoint to avoid overlap
                if xt(i) >= x(j) && xt(i) < x(j+1)
                    % Calculate offset from left endpoint of segment
                    xh = xt(i) - x(j);
                    
                    % Evaluate cubic polynomial: a + b*xh + c*xh^2 + d*xh^3
                    g(i) = C(4*j-3) + C(4*j-2)*xh + C(4*j-1)*xh^2 + C(4*j)*xh^3;
                    break;
                end
            end
        end
    end

end

function g = evaluate_spline(x_data, coeffs, xt)
% evaluate_spline - Evaluate cubic spline at given points
%
% Inputs:
% x_data - original knot points
% coeffs - spline coefficients from cubic_spline function
% xt     - evaluation points
%
% Output:
% g      - spline values at xt points

    n = length(x_data) - 1;
    m = length(xt);
    g = zeros(size(xt));
    
    for i = 1:m
        % Find which segment contains the evaluation point
        if xt(i) <= x_data(1)
            j = 1;  % Extrapolate using first segment
        elseif xt(i) >= x_data(end)
            j = n;  % Extrapolate using last segment
        else
            j = find(x_data <= xt(i), 1, 'last');
            if j > n, j = n; end
        end
        
        % Calculate offset from left endpoint of segment
        xh = xt(i) - x_data(j);
        coeff_idx = 4*j - 3;
        
        % Evaluate cubic polynomial: S_j(x) = a + b*xh + c*xh^2 + d*xh^3
        % Using Horner's method for numerical stability
        g(i) = coeffs(coeff_idx) + xh * (coeffs(coeff_idx+1) + ...
               xh * (coeffs(coeff_idx+2) + xh * coeffs(coeff_idx+3)));
    end
end

function bc_info = analyze_boundary_conditions(x_data, y_data)
% analyze_boundary_conditions - Analyze appropriate boundary conditions
% for gamma function data and provide justification
%
% Inputs:
% x_data, y_data - data points
%
% Output:
% bc_info - structure with boundary condition analysis

    bc_info = struct();
    
    % Natural spline analysis
    bc_info.natural.description = 'Natural spline: S''''(x1) = S''''(x5) = 0';
    bc_info.natural.justification = ['No derivative information available. ', ...
        'Natural splines minimize curvature and provide smooth interpolation.'];
    bc_info.natural.bc = [0, 0, 2];
    
    % Estimate derivatives for clamped spline using finite differences
    % Left endpoint derivative (forward difference)
    h1 = x_data(2) - x_data(1);
    h2 = x_data(3) - x_data(2);
    left_deriv = (-3*y_data(1) + 4*y_data(2) - y_data(3)) / (2*h1);
    
    % Right endpoint derivative (backward difference)
    hn_1 = x_data(end-1) - x_data(end-2);
    hn = x_data(end) - x_data(end-1);
    right_deriv = (y_data(end-2) - 4*y_data(end-1) + 3*y_data(end)) / (2*hn);
    
    bc_info.clamped.description = 'Clamped spline with estimated derivatives';
    bc_info.clamped.justification = ['Estimated derivatives using finite differences. ', ...
        'Provides better end behavior when derivative estimates are reasonable.'];
    bc_info.clamped.left_deriv = left_deriv;
    bc_info.clamped.right_deriv = right_deriv;
    bc_info.clamped.bc = [left_deriv, right_deriv, 1];
    
    % Theoretical gamma function derivatives (for comparison)
    bc_info.theoretical.description = 'Theoretical gamma function derivatives';
    bc_info.theoretical.justification = ['Uses known properties: Gamma''(1) ≈ -0.5772 ', ...
        '(Euler constant), Gamma''(5) estimated from gamma function behavior.'];
    bc_info.theoretical.bc = [-0.5772, 52.34, 1];
end