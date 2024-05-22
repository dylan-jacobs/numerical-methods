% Composite Gauss-Legendre 

function [approximation] = CompositeGaussLegendre(a, b, N, M, f)
    % interval [a, b]
    % N intervals = N quadratures
    % M nodes, weights
    % function f

    deltaX = (b - a) / N;
    approximation = 0;
    for i = 1:N
        sub_a = a + (deltaX * (i - 1));
        sub_b = a + (deltaX * (i));
        [weights, nodes] = GetWeightsAndNodes(M);
        sub_approx = GaussLegendre(sub_a, sub_b, weights, nodes, f);
        approximation = approximation + sub_approx;
    end
end