function [M, T, S, G] = Nint(a, b, f, n)
    M = 0;
    T = 0;
    S = 0;
    G = 0;
    h = (b-a)/n;
    x = a:h:b;
    
    for i = 1:n
        % Midpoint Rule - Fixed: removed extra parenthesis causing division issue
        M = M + h * feval(f, (x(i) + x(i+1))/2);
        
        % Trapezoidal Rule - Fixed: moved h/2 multiplication outside the sum
        T = T + (h/2) * (feval(f, x(i)) + feval(f, x(i+1)));
        
        % Simpson's Rule - Fixed: check if i is odd/even, not mod(1,2)
        if mod(i, 2) == 1  % odd index
            if i == 1
                S = S + (h/3) * feval(f, x(i));
            else
                S = S + (2*h/3) * feval(f, x(i));
            end
        else  % even index
            S = S + (4*h/3) * feval(f, x(i));
        end
    end
    
    % Add the final point for Simpson's rule
    S = S + (h/3) * feval(f, x(n+1));
end