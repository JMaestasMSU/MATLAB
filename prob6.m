A = [0 1 0 0; 1 0 3 0; -0.5 0 -0.2 1; -0.5 -0.3 1 0];

% Step 1 (we already did this by hand)
P1 = [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
A1 = P1 * A

M1 = [1 0 0 0; 0 1 0 0; 0.5 0 1 0; 0.5 0 0 1];
A2 = M1 * A1

% Step 2 - Find M2
fprintf('Looking at positions (2,2) through (4,2): [%.1f, %.1f, %.1f]\n', A2(2,2), A2(3,2), A2(4,2));
fprintf('Largest magnitude is %.1f in row 2\n', A2(2,2));
fprintf('No pivot needed.\n\n');

P2 = eye(4);

% Compute multiplier for row 4 (only nonzero below pivot)
m42 = A2(4,2) / A2(2,2);
fprintf('Multiplier m42 = %.1f / %.1f = %.1f\n', A2(4,2), A2(2,2), m42);

M2 = eye(4);
M2(4,2) = -m42;  % Store negative of multiplier
fprintf('\nM^(2) = \n');
disp(M2);

A3 = M2 * A2;
fprintf('After M^(2)A^(2):\n');
disp(A3);

% Step 3 - Find M3
fprintf('Looking at positions (3,3) through (4,3): [%.1f, %.1f]\n', A3(3,3), A3(4,3));
[~, pivot_row] = max(abs(A3(3:4, 3)));
pivot_row = pivot_row + 2;  % Adjust for actual row number
fprintf('Largest magnitude is %.1f in row %d\n', abs(A3(pivot_row,3)), pivot_row);
fprintf('Pivot needed. Exchange rows 3 and 4.\n\n');

P3 = eye(4);
P3([3 4], :) = P3([4 3], :);  % Swap rows 3 and 4
fprintf('P^(3) = \n');
disp(P3);

A4 = P3 * A3;
fprintf('After P^(3)A^(3):\n');
disp(A4);

% Compute multiplier for row 4
m43 = A4(4,3) / A4(3,3);
fprintf('Multiplier m43 = %.2f / %.1f = %.2f\n', A4(4,3), A4(3,3), m43);

M3 = eye(4);
M3(4,3) = -m43;  % Store negative of multiplier
fprintf('\nM^(3) = \n');
disp(M3);

U = M3 * A4;
fprintf('Final U = \n');
disp(U);

% Compute M-tilde matrices
M1_tilde = M1 * P2 * P3;
fprintf('M_tilde^(1) = M^(1)P^(2)P^(3) = \n');
disp(M1_tilde);

M2_tilde = M2 * P3;
fprintf('M_tilde^(2) = M^(2)P^(3) = \n');
disp(M2_tilde);

M3_tilde = M3;
fprintf('M_tilde^(3) = M^(3) = \n');
disp(M3_tilde);

% Compute L
fprintf('Computing L\n');
L = inv(M1_tilde) * inv(M2_tilde) * inv(M3_tilde);
fprintf('L = [M_tilde^(1)]^(-1) [M_tilde^(2)]^(-1) [M_tilde^(3)]^(-1) = \n');
disp(L);

% Compute overall P
P = P3 * P2 * P1;
fprintf('P = P^(3)P^(2)P^(1) = \n');
disp(P);

% Verify
fprintf('||PA - LU|| = %.2e\n', norm(P*A - L*U, 'fro'));