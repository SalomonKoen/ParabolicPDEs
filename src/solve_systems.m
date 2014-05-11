function [UU, VV] = solve_systems(D_u, D_v, a, b, L, a_x, b_x, M, N, k, U_c, V_c)
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
%   U_c   := Vector for initial prey distribution.
%   G   := Vector for initial predator distribution.
%
% Output:
%
%   UU := Prey over time. (time steps as columns)
%   VV := Predators over time. (time steps as columns)
%

%%% Preallocations.

% Solution matrices "UU" and "VV".
UU = zeros(N,M);
VV = zeros(N,M);

% Current and next time step solution vectors for prey. Initial values
% are random. (c = current, n = next)
U_n = zeros(N,1);
UU(:,1) = U_c;

% Similarly for predators.
V_n = zeros(N,1);
VV(:,1) = V_c;

%%%

% Setup waitbar.
wb = waitbar(0, sprintf('%3.2f%%', 0), 'Name', 'Working...',...
    'CreateCancelBtn',...
    'setappdata(gcbf, ''canceling'',1)');

setappdata(wb, 'canceling', 0);

% Solve system, filling "UU" and "VV" columns from left to right.
for m = 2:M
    % Escape computation when cancel button is pressed.
    % The partial result is returned for plotting.
    if getappdata(wb, 'canceling')
        break;
    end
    
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
    
    % Update waitbar.
    waitbar(m / M, wb, sprintf('%3.2f%%', (m / M) * 100));
end

% Remove the waitbar.
delete(wb);

end

