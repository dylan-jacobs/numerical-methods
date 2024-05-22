% Gaussian Quadatures
% Numerical methods for solving definite integrals

%% Convergence Test 1
clc;
clear variables;
close all;

a = 0; b = 2; N = 3; M = 4;
f = @(x) sin(x);
[weights, nodes] = GetWeightsAndNodes(N);

approx = GaussLegendre(a, b, weights, nodes, f)
exact = -cos(2) + cos(0)

comp_approx = CompositeGaussLegendre(a, b, N, M, f)
comp_exact =-cos(2) + cos(0)

%% Convergence Test 2

a = 0; b = 2; 
NVals = [4; 8; 16; 32; 64; 128; 256]; 
MVals = [2, 3];
f = @(x) sin(x);
exact_integral = -cos(b) + cos(a);

NVals = [4; 8; 16; 32; 64; 128; 256];

approximations = zeros(2, numel(NVals), 1);
for i = 1:numel(NVals)
    N = NVals(i);
    approximations(1, i) = CompositeGaussLegendre(a, b, N, MVals(1), f);
    approximations(2, i) = CompositeGaussLegendre(a, b, N, MVals(2), f);
end
order2 = CreateErrorPlot(a, b, exact_integral, approximations(1, :), NVals, 4, 'sin(x)')'
order6 = CreateErrorPlot(a, b, exact_integral, approximations(2, :), NVals, 6, 'sin(x)')'



%% Convergence Test 3
a = -1; b = 2; N = 3; 
MVals = [2, 3];
f = @(x) exp(-x .^ 2);
exact_integral = (sqrt(pi)/2)*(erf(b)-erf(a));

NVals = [4; 8; 16; 32; 64; 128; 256];

approximations = zeros(2, numel(NVals), 1);
for i = 1:numel(NVals)
    N = NVals(i);
    approximations(1, i) = CompositeGaussLegendre(a, b, N, MVals(1), f);
    approximations(2, i) = CompositeGaussLegendre(a, b, N, MVals(2), f);
end
order2 = CreateErrorPlot(a, b, exact_integral, approximations(1, :), NVals, 4, 'exp(-x .^ 2)')'
order6 = CreateErrorPlot(a, b, exact_integral, approximations(2, :), NVals, 6, 'exp(-x .^ 2)')'

































