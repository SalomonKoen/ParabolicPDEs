function [prey, pred] = init_vectors(d, r, m, a, b, x_a, x_b, N)
%INIT_VECTORS Build initial predator and prey vectors.
%
% Input:
%   d := Prey population density.
%   r := Prey randomness factor for density.
%
%   m := Magnitude of predator invasion front.
%   a := Left end point of predator invasion front.
%   b := Right end point of predator invasion front.
%
%   x_a := Left end point of x domain.
%   x_b := Right end point of x domain.
%
%   N := Length of output vectors.
%
% Output:
%   prey := Resultant prey vector for t=0.
%   pred := Resultant predator vector for t=0.

% Allocate output.
prey = d + (r * rand(N,1));
pred = linspace(x_a, x_b, N)';

% Build stepped invasion front.
for i = 1:N
    if (pred(i) >= a) && (pred(i) <= b)
        pred(i) = m;
    else
        pred(i) = 0;
    end
end

end

