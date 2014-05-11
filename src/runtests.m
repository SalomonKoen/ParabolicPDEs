function [] = runtests()
%RUNTESTS Execute set of tests for Crank-Nicolson model.
%   Each test should be written in a seperate subfunction
%   in this file.

clear; clc; % Clean workspace first.

test_solve_systems();

end

function [success] = test_solve_systems()

% Parameters are very sensitive.
D_u = 0.09;
D_v = 0.2;
a   = 3;
b   = 1;
L   = 1;
a_x = 0;
b_x = 3;
M   = 1000;
N   = 35;
k   = 0.009;
F   = @(x) 1+x*0.1; % Watch out for negative values here.
G   = @(x) 0.9-0.4*x; % And here. "a_x" and "b_x" should be chosen carefully.

[UU, VV] = solve_systems(D_u, D_v, a, b, L, a_x, b_x, M, N, k, F, G);

rotate3d on
hold on
surf(UU, 'EdgeColor', 'none')
surf(VV, 'EdgeColor', 'none')
hold off

success = true;

end