% Define the matrix
A = [-2 1 0 0; 1 -3 -1 0; 0 -1 1 1; 0 0 1 3];

fprintf('========== ORIGINAL MATRIX ==========\n');
fprintf('A^(1) =\n');
disp(A)

%% ========== ITERATION 1 WITHOUT SHIFTING ==========

fprintf('\n');
fprintf('================================================================================\n');
fprintf('                    ITERATION 1 - WITHOUT SHIFTING\n');
fprintf('================================================================================\n\n');

% Step 1: First Givens rotation to eliminate (2,1)
fprintf('---------- STEP 1: Build P_2 to eliminate (2,1) ----------\n\n');
fprintf('Column 1 of A^(1): [%.6f; %.6f; %.6f; %.6f]\n\n', A(1,1), A(2,1), A(3,1), A(4,1));

a11 = A(1,1);
a21 = A(2,1);
r1 = sqrt(a11^2 + a21^2);
c1 = a11 / r1;
s1 = a21 / r1;

fprintf('r = sqrt(%.6f^2 + %.6f^2) = %.6f\n', a11, a21, r1);
fprintf('c = %.6f / %.6f = %.6f\n', a11, r1, c1);
fprintf('s = %.6f / %.6f = %.6f\n\n', a21, r1, s1);

P2 = eye(4);
P2(1,1) = c1;
P2(1,2) = -s1;
P2(2,1) = s1;
P2(2,2) = c1;

fprintf('P_2 =\n');
for i = 1:4
    fprintf('[%.6f  %.6f  %.6f  %.6f]\n', P2(i,:));
end
fprintf('\n');

P2_A = P2 * A;
fprintf('P_2 * A^(1) =\n');
for i = 1:4
    fprintf('[%.6f  %.6f  %.6f  %.6f]\n', P2_A(i,:));
end
fprintf('\n');

% Step 2: Second Givens rotation to eliminate (3,2)
fprintf('---------- STEP 2: Build P_3 to eliminate (3,2) ----------\n\n');
fprintf('Column 2 (rows 2-3) of P_2*A^(1): [%.6f; %.6f]\n\n', P2_A(2,2), P2_A(3,2));

a22 = P2_A(2,2);
a32 = P2_A(3,2);
r2 = sqrt(a22^2 + a32^2);
c2 = a22 / r2;
s2 = a32 / r2;

fprintf('r = sqrt(%.6f^2 + %.6f^2) = %.6f\n', a22, a32, r2);
fprintf('c = %.6f / %.6f = %.6f\n', a22, r2, c2);
fprintf('s = %.6f / %.6f = %.6f\n\n', a32, r2, s2);

P3 = eye(4);
P3(2,2) = c2;
P3(2,3) = -s2;
P3(3,2) = s2;
P3(3,3) = c2;

fprintf('P_3 =\n');
for i = 1:4
    fprintf('[%.6f  %.6f  %.6f  %.6f]\n', P3(i,:));
end
fprintf('\n');

P3_P2_A = P3 * P2_A;
fprintf('P_3 * P_2 * A^(1) =\n');
for i = 1:4
    fprintf('[%.6f  %.6f  %.6f  %.6f]\n', P3_P2_A(i,:));
end
fprintf('\n');

% Step 3: Third Givens rotation to eliminate (4,3)
fprintf('---------- STEP 3: Build P_4 to eliminate (4,3) ----------\n\n');
fprintf('Column 3 (rows 3-4) of P_3*P_2*A^(1): [%.6f; %.6f]\n\n', P3_P2_A(3,3), P3_P2_A(4,3));

a33 = P3_P2_A(3,3);
a43 = P3_P2_A(4,3);
r3 = sqrt(a33^2 + a43^2);
c3 = a33 / r3;
s3 = a43 / r3;

fprintf('r = sqrt(%.6f^2 + %.6f^2) = %.6f\n', a33, a43, r3);
fprintf('c = %.6f / %.6f = %.6f\n', a33, r3, c3);
fprintf('s = %.6f / %.6f = %.6f\n\n', a43, r3, s3);

P4 = eye(4);
P4(3,3) = c3;
P4(3,4) = -s3;
P4(4,3) = s3;
P4(4,4) = c3;

fprintf('P_4 =\n');
for i = 1:4
    fprintf('[%.6f  %.6f  %.6f  %.6f]\n', P4(i,:));
end
fprintf('\n');

% Step 4: Compute R^(1)
fprintf('---------- STEP 4: Compute R^(1) = P_4 * P_3 * P_2 * A^(1) ----------\n\n');

R1 = P4 * P3_P2_A;
fprintf('R^(1) =\n');
for i = 1:4
    fprintf('[%.6f  %.6f  %.6f  %.6f]\n', R1(i,:));
end
fprintf('\n');

% Step 5: Compute Q^(1)
fprintf('---------- STEP 5: Compute Q^(1) = P_2^T * P_3^T * P_4^T ----------\n\n');

Q1 = P2' * P3' * P4';
fprintf('Q^(1) =\n');
for i = 1:4
    fprintf('[%.6f  %.6f  %.6f  %.6f]\n', Q1(i,:));
end
fprintf('\n');

% Verify
fprintf('---------- VERIFICATION ----------\n\n');
fprintf('Q^(1) * R^(1) = A^(1)?\n');
QR = Q1 * R1;
for i = 1:4
    fprintf('[%.6f  %.6f  %.6f  %.6f]\n', QR(i,:));
end
fprintf('\n');

% Step 6: Compute A^(2)
fprintf('---------- STEP 6: Compute A^(2) = R^(1) * Q^(1) ----------\n\n');

A2 = R1 * Q1;
fprintf('A^(2) =\n');
for i = 1:4
    fprintf('[%.6f  %.6f  %.6f  %.6f]\n', A2(i,:));
end
fprintf('\n');

%% ========== ITERATION 2 WITHOUT SHIFTING ==========

fprintf('\n');
fprintf('================================================================================\n');
fprintf('                    ITERATION 2 - WITHOUT SHIFTING\n');
fprintf('================================================================================\n\n');

fprintf('Starting with A^(2) =\n');
for i = 1:4
    fprintf('[%.6f  %.6f  %.6f  %.6f]\n', A2(i,:));
end
fprintf('\n');

% Step 1: First Givens rotation to eliminate (2,1)
fprintf('---------- STEP 1: Build P_2 to eliminate (2,1) ----------\n\n');

a11_2 = A2(1,1);
a21_2 = A2(2,1);
r1_2 = sqrt(a11_2^2 + a21_2^2);
c1_2 = a11_2 / r1_2;
s1_2 = a21_2 / r1_2;

fprintf('r = sqrt(%.6f^2 + %.6f^2) = %.6f\n', a11_2, a21_2, r1_2);
fprintf('c = %.6f / %.6f = %.6f\n', a11_2, r1_2, c1_2);
fprintf('s = %.6f / %.6f = %.6f\n\n', a21_2, r1_2, s1_2);

P2_iter2 = eye(4);
P2_iter2(1,1) = c1_2;
P2_iter2(1,2) = -s1_2;
P2_iter2(2,1) = s1_2;
P2_iter2(2,2) = c1_2;

P2_A2 = P2_iter2 * A2;

% Step 2: Second Givens rotation to eliminate (3,2)
fprintf('---------- STEP 2: Build P_3 to eliminate (3,2) ----------\n\n');

a22_2 = P2_A2(2,2);
a32_2 = P2_A2(3,2);
r2_2 = sqrt(a22_2^2 + a32_2^2);
c2_2 = a22_2 / r2_2;
s2_2 = a32_2 / r2_2;

fprintf('r = sqrt(%.6f^2 + %.6f^2) = %.6f\n', a22_2, a32_2, r2_2);
fprintf('c = %.6f / %.6f = %.6f\n', a22_2, r2_2, c2_2);
fprintf('s = %.6f / %.6f = %.6f\n\n', a32_2, r2_2, s2_2);

P3_iter2 = eye(4);
P3_iter2(2,2) = c2_2;
P3_iter2(2,3) = -s2_2;
P3_iter2(3,2) = s2_2;
P3_iter2(3,3) = c2_2;

P3_P2_A2 = P3_iter2 * P2_A2;

% Step 3: Third Givens rotation to eliminate (4,3)
fprintf('---------- STEP 3: Build P_4 to eliminate (4,3) ----------\n\n');

a33_2 = P3_P2_A2(3,3);
a43_2 = P3_P2_A2(4,3);
r3_2 = sqrt(a33_2^2 + a43_2^2);
c3_2 = a33_2 / r3_2;
s3_2 = a43_2 / r3_2;

fprintf('r = sqrt(%.6f^2 + %.6f^2) = %.6f\n', a33_2, a43_2, r3_2);
fprintf('c = %.6f / %.6f = %.6f\n', a33_2, r3_2, c3_2);
fprintf('s = %.6f / %.6f = %.6f\n\n', a43_2, r3_2, s3_2);

P4_iter2 = eye(4);
P4_iter2(3,3) = c3_2;
P4_iter2(3,4) = -s3_2;
P4_iter2(4,3) = s3_2;
P4_iter2(4,4) = c3_2;

% Compute R^(2) and Q^(2)
R2 = P4_iter2 * P3_P2_A2;
Q2 = P2_iter2' * P3_iter2' * P4_iter2';

fprintf('---------- STEP 4: Compute A^(3) = R^(2) * Q^(2) ----------\n\n');

A3 = R2 * Q2;
fprintf('A^(3) =\n');
for i = 1:4
    fprintf('[%.6f  %.6f  %.6f  %.6f]\n', A3(i,:));
end
fprintf('\n');

%% ========== ITERATION 1 WITH SHIFTING ==========

fprintf('\n');
fprintf('================================================================================\n');
fprintf('                    ITERATION 1 - WITH SHIFTING\n');
fprintf('================================================================================\n\n');

shift1 = A(4,4);
fprintf('Choose shift: sigma_1 = A^(1)(4,4) = %.6f\n\n', shift1);

A_shifted1 = A - shift1 * eye(4);
fprintf('A^(1) - sigma_1*I =\n');
for i = 1:4
    fprintf('[%.6f  %.6f  %.6f  %.6f]\n', A_shifted1(i,:));
end
fprintf('\n');

% Build Givens rotations
a11_s = A_shifted1(1,1);
a21_s = A_shifted1(2,1);
r1_s = sqrt(a11_s^2 + a21_s^2);
c1_s = a11_s / r1_s;
s1_s = a21_s / r1_s;

P2_shift1 = eye(4);
P2_shift1(1,1) = c1_s;
P2_shift1(1,2) = -s1_s;
P2_shift1(2,1) = s1_s;
P2_shift1(2,2) = c1_s;

P2_A_shift1 = P2_shift1 * A_shifted1;

a22_s = P2_A_shift1(2,2);
a32_s = P2_A_shift1(3,2);
r2_s = sqrt(a22_s^2 + a32_s^2);
c2_s = a22_s / r2_s;
s2_s = a32_s / r2_s;

P3_shift1 = eye(4);
P3_shift1(2,2) = c2_s;
P3_shift1(2,3) = -s2_s;
P3_shift1(3,2) = s2_s;
P3_shift1(3,3) = c2_s;

P3_P2_A_shift1 = P3_shift1 * P2_A_shift1;

a33_s = P3_P2_A_shift1(3,3);
a43_s = P3_P2_A_shift1(4,3);
r3_s = sqrt(a33_s^2 + a43_s^2);
c3_s = a33_s / r3_s;
s3_s = a43_s / r3_s;

P4_shift1 = eye(4);
P4_shift1(3,3) = c3_s;
P4_shift1(3,4) = -s3_s;
P4_shift1(4,3) = s3_s;
P4_shift1(4,4) = c3_s;

R_shift1 = P4_shift1 * P3_P2_A_shift1;
Q_shift1 = P2_shift1' * P3_shift1' * P4_shift1';

A2_shift1 = R_shift1 * Q_shift1 + shift1 * eye(4);
fprintf('A^(2) (with shift) =\n');
for i = 1:4
    fprintf('[%.6f  %.6f  %.6f  %.6f]\n', A2_shift1(i,:));
end
fprintf('\n');

%% ========== ITERATION 2 WITH SHIFTING ==========

fprintf('\n');
fprintf('================================================================================\n');
fprintf('                    ITERATION 2 - WITH SHIFTING\n');
fprintf('================================================================================\n\n');

shift2 = A2_shift1(4,4);
fprintf('Choose shift: sigma_2 = A^(2)(4,4) = %.6f\n\n', shift2);

A_shifted2 = A2_shift1 - shift2 * eye(4);

% Build Givens rotations (same process as before)
a11_s2 = A_shifted2(1,1);
a21_s2 = A_shifted2(2,1);
r1_s2 = sqrt(a11_s2^2 + a21_s2^2);
c1_s2 = a11_s2 / r1_s2;
s1_s2 = a21_s2 / r1_s2;

P2_shift2 = eye(4);
P2_shift2(1,1) = c1_s2;
P2_shift2(1,2) = -s1_s2;
P2_shift2(2,1) = s1_s2;
P2_shift2(2,2) = c1_s2;

P2_A_shift2 = P2_shift2 * A_shifted2;

a22_s2 = P2_A_shift2(2,2);
a32_s2 = P2_A_shift2(3,2);
r2_s2 = sqrt(a22_s2^2 + a32_s2^2);
c2_s2 = a22_s2 / r2_s2;
s2_s2 = a32_s2 / r2_s2;

P3_shift2 = eye(4);
P3_shift2(2,2) = c2_s2;
P3_shift2(2,3) = -s2_s2;
P3_shift2(3,2) = s2_s2;
P3_shift2(3,3) = c2_s2;

P3_P2_A_shift2 = P3_shift2 * P2_A_shift2;

a33_s2 = P3_P2_A_shift2(3,3);
a43_s2 = P3_P2_A_shift2(4,3);
r3_s2 = sqrt(a33_s2^2 + a43_s2^2);
c3_s2 = a33_s2 / r3_s2;
s3_s2 = a43_s2 / r3_s2;

P4_shift2 = eye(4);
P4_shift2(3,3) = c3_s2;
P4_shift2(3,4) = -s3_s2;
P4_shift2(4,3) = s3_s2;
P4_shift2(4,4) = c3_s2;

R_shift2 = P4_shift2 * P3_P2_A_shift2;
Q_shift2 = P2_shift2' * P3_shift2' * P4_shift2';

A3_shift2 = R_shift2 * Q_shift2 + shift2 * eye(4);
fprintf('A^(3) (with shift) =\n');
for i = 1:4
    fprintf('[%.6f  %.6f  %.6f  %.6f]\n', A3_shift2(i,:));
end
fprintf('\n');

%% SUMMARY

fprintf('\n');
fprintf('================================================================================\n');
fprintf('                              SUMMARY\n');
fprintf('================================================================================\n\n');

fprintf('Actual eigenvalues:\n');
evals = sort(eig(A), 'descend');
fprintf('  %.6f\n', evals);

fprintf('\nDiagonal of A^(2) without shifting:\n');
diag_A2 = sort(diag(A2), 'descend');
fprintf('  %.6f\n', diag_A2);

fprintf('\nDiagonal of A^(3) without shifting:\n');
diag_A3 = sort(diag(A3), 'descend');
fprintf('  %.6f\n', diag_A3);

fprintf('\nDiagonal of A^(2) with shifting:\n');
diag_A2_shift = sort(diag(A2_shift1), 'descend');
fprintf('  %.6f\n', diag_A2_shift);

fprintf('\nDiagonal of A^(3) with shifting:\n');
diag_A3_shift = sort(diag(A3_shift2), 'descend');
fprintf('  %.6f\n', diag_A3_shift);