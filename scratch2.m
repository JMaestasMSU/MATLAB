% Set up
x0 = 0;
exact_value = exp(0);  % f'''(0) = e^0 = 1

% h values to test
h_values = 10.^(-1:-1:-9);  % [10^-1, 10^-2, ..., 10^-9]

% Storage
approximations = zeros(size(h_values));
errors = zeros(size(h_values));

% Loop through each h
for i = 1:length(h_values)
    h = h_values(i);
    
    % Evaluate function at the four points
    f_plus2h = exp(x0 + 2*h);
    f_plush = exp(x0 + h);
    f_minush = exp(x0 - h);
    f_minus2h = exp(x0 - 2*h);
    
    % Apply our formula
    f_triple_prime = (f_plus2h - 2*f_plush + 2*f_minush - f_minus2h) / (2*h^3);
    
    % Store results
    approximations(i) = f_triple_prime;
    errors(i) = abs(f_triple_prime - exact_value);
end