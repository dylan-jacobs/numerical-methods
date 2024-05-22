% Second order centered difference
% Discretize u''(x) on [a, b]

clc;
clear variables;
close all;

a = 0; b = 2 * pi; % interval [a, b]
Nvals = [4, 8, 16, 32, 64, 128]; % N cells, N+1 points

L_inf_vals = zeros(numel(Nvals), 1);
L1_vals = zeros(numel(Nvals), 1);
L2_vals = zeros(numel(Nvals), 1);

for k = 1:numel(Nvals)
    N = Nvals(k);
    xvals = linspace(a, b, N + 1)';
    uvals = exp(sin(xvals)); 
    dx = xvals(2) - xvals(1);
    
    uxx_exact = exp(sin(xvals)) .* (cos(xvals) .^ 2 - sin(xvals));
    
    % construct differentiation matrix for u''(x)
    Dxx = (gallery('tridiag', (N + 1), 1, -2, 1));
    Dxx(1, end - 1) = 1; Dxx(end, 2) = 1; % boundary conditions
    Dxx = Dxx / (dx ^ 2);
    
    uxx_approx = Dxx * uvals;
    L_inf_vals(k) = max(abs(uxx_exact - uxx_approx));
    L1_vals(k) = dx * sum(abs(uxx_exact - uxx_approx));
    L2_vals(k) = sqrt(sum((uxx_exact - uxx_approx) .^ 2) * dx);
end

log2(L_inf_vals(1:end-1) ./ L_inf_vals(2:end));
log2(L1_vals(1:end-1) ./ L1_vals(2:end));
log2(L2_vals(1:end-1) ./ L2_vals(2:end));

%% HOMEWORK
% Use second-order difference formula 
% I = [0, pi], u(0) = 0, u(pi) = 0
% u = sin(x) / (1 + (x-(pi/2)^2))
% only use formula on interior nodes because we know solution at 0, pi
% 


% Find u(x) to test this formula on [0, pi] where u(0) = 1, u(pi) = 0

































