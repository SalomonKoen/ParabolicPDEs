function [UU, VV] = solve_systems(D_u, D_v, a, b, L, a_x, b_x, M, N, k, F, G)
%SOLVE_SYSTEMS Solve, via Crank-Nicolson, the coupled diffusion PDEs.
%
% Input:
%
%   D_u := Prey diffusion.
%   D_v := Predator diffusion.
%
%   a   := Model parameter.
%   b   := Model parameter.
%   L   := Model parameter (\lambda).
%
%   a_x := Left boundary.
%   b_x := Right boundary.
%
%   M   := Number of time steps for temporal domain.
%   N   := Number of points in spatial domain to solve for.
%
%   k   := Step length for temporal domain.
%
%   F   := Function for initial prey distribution.
%   G   := Function for initial predator distribution.
%
%     Functions "F" and "G" take single parameter "x". Should use
%     "@(x) rand(length(x),1)" for random initial values.
%
% Output:
%
%   UU := Prey over time. (time steps as columns)
%   VV := Predators over time. (time steps as columns)
%

%%% Preallocations.

X = linspace(a_x, b_x, N)'; % Must transpose to column vector.

% Solution matrices "UU" and "VV".
UU = zeros(N,M);
VV = zeros(N,M);

% Current and next time step solution vectors for prey. Initial values
% are random. (c = current, n = next)
U_c = feval(F, X);
U_n = zeros(N,1);

UU(:,1) = U_c;

% Similarly for predators.
V_c = feval(G, X);
V_n = zeros(N,1);

VV(:,1) = V_c;

%%%

% Solve system, filling "UU" and "VV" columns from left to right.
for m = 2:M
    % Construct matrices for current time step.
    [U_L, U_R, V_L, V_R] = build_matrices(D_u, D_v, a, b, L, a_x, b_x, N, k, U_c, V_c);
    
    % Solve both systems.
    U_n = U_L \ (U_R * U_c);
    V_n = V_L \ (V_R * V_c);
    
    % Store solutions to current time step in output matrices.
    UU(:,m) = U_n;
    VV(:,m) = V_n;
    
    % Replace current with next vectors.
    U_c = U_n;
    V_c = V_n;
end

end

