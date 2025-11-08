x = (0:0.25:1)';  % column vector
y = 0:0.25:1;     % row vector
D = x - y;        % outer subtraction gives 5×5 matrix
F = exp(D);