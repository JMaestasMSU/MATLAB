% GAUSS  Compute nodes x (Legendre points) and weights w for Gauss quadrature
%        using the Golub-Welsch algorithm.
%        Code from book Spectral Methods by L. Trefethen
%
%   [x, w] = gauss(N)
%
%   Inputs:
%     N - Number of quadrature points (integer, N >= 1)
%
%   Outputs:
%     x - Column vector of N Gauss-Legendre nodes (roots of Legendre polynomial)
%     w - Column vector of N corresponding quadrature weights
%
%   Example:
%     [x, w] = gauss(4)
%     % Returns 4-point Gauss-Legendre nodes and weights on [-1, 1]
%
%   This function constructs the symmetric tridiagonal Jacobi matrix
%   whose eigenvalues are the Legendre nodes and whose eigenvectors
%   yield the weights.

function [x, w] = gauss(N)
    % Compute the off-diagonal elements (beta) for the Jacobi matrix
    beta = .5 ./ sqrt(1 - (2*(1:N-1)).^(-2));
    
    % Construct the symmetric tridiagonal Jacobi matrix T
    T = diag(beta, 1) + diag(beta, -1);
    
    % Compute eigenvalues (nodes) and eigenvectors of T
    [V, D] = eig(T);
    
    % Extract and sort the eigenvalues (nodes)
    x = diag(D);
    [x, i] = sort(x);
    
    % Compute the weights from the first row of eigenvectors
    w = 2 * V(1, i).^2;
end