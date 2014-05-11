function [U_L, U_R, V_L, V_R] = build_matrices(D_u, D_v, a, b, L, a_x, b_x, N, k, U, V)
%BUILD_MATRICES Construct matrices to solve Crank-Nicolson model.
%
% Input:
%
%   D_u := Prey diffusion.
%   D_u := Predator diffusion.
%
%   a   := Model parameter.
%   b   := Model parameter.
%   L   := Model parameter (\lambda).
%
%   a_x := Left boundary.
%   b_x := Right boundary.
%   N   := Number of points in x domain.
%
%   k   := Step length for for t domain.
%
%   U   := Prey vector on current time step.
%   V   := Predator vector on current time step.
%
% Output:
%
%   U_L := Sparse matrix for prey system LHS.
%   U_R := Sparse matrix for prey system RHS.
%   V_L := Sparse matrix for predator system LHS.
%   V_R := Sparse matrix for predator system RHS.
%

%%% Precompute scalar values.
h = (b_x - a_x) / (N - 1);

r = (D_u * k) / h^2;
s = (D_v * k) / h^2;
%%%

%%% Main diagonal expressions. (vectorised)

% Prey system.
L_u = 2 + (2 * r) - k + (k .* U) + (a * k .* V) ./ (1 + L .* U);
R_u = 2 - (2 * r) + k - (k .* U) - (a * k .* V) ./ (1 + L .* U);

R = r * ones(N, 1);

U_L = spdiags([-R L_u -R], -1:1, N, N);
U_R = spdiags([ R R_u  R], -1:1, N, N);

% Predator system.
L_v = 2 + (2 * s) + (k / (a * b)) - ((a * k .* U) ./ (b * (1 + L .* U)));
R_v = 2 - (2 * s) - (k / (a * b)) + ((a * k .* U) ./ (b * (1 + L .* U)));

S = s * ones(N, 1);

V_L = spdiags([-S L_v -S], -1:1, N, N);
V_R = spdiags([ S R_v  S], -1:1, N, N);

%%%

%%% Neumann conditions at endpoints "a_x" and "b_x". (double up values)

% Prey system.
U_L(1,2)   = -2*r;
U_L(N,N-1) = -2*r;

U_R(1,2)   =  2*r;
U_R(N,N-1) =  2*r;

% Predator system.
V_L(1,2)   = -2*s;
V_L(N,N-1) = -2*s;

V_R(1,2)   =  2*s;
V_R(N,N-1) =  2*s;

%%%

end
