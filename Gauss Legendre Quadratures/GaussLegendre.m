% Single Gauss-Legendre Quadrature approximation
% We do a linear mapping from [-1, 1] to [a, b]
% Use greek letter Xi to correspond to integral from [a, b] --> f(Xi)
% x = -1 + (2/(b-a))(Xi - a)
% dx = (2/(b-a))dXi

function [approximation] = GaussLegendre(a, b, weights, nodes, f)
    % interval [a, b]
    % N intervals, N+1 nodes
    % f = function

    approximation = 0;
    deltaXi = (b-a) / 2;
    for i = 1:numel(nodes)
        node = nodes(i);
        Xi = ((node + 1) * (b-a) / 2) + a;
        Xi = node * (b-a) / 2;
        approximation = approximation + (weights(i) * f(Xi + (b+a)/2));
    end
    approximation = deltaXi * approximation;
end
