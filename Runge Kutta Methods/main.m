% main script to test all Runge-Kutta approximations
clc;
clear variables;
close all;
 
% function f
f = @(u, t) 2*u*t;
f_exact = @(t) exp(t^2);

% initial condition
u0 = 1;

interval = [0, 0.5]; % interval over which to approximate function f
Nvals = [4, 8, 16, 32, 64, 128, 256];

% compute exact value of f at final time
exact = f_exact(interval(2));

% calculate errors
expRK4_approx = zeros(numel(Nvals), 1);
sspRK3_approx = zeros(numel(Nvals), 1);
sspRK2_approx = zeros(numel(Nvals), 1);

% compute errors for each N value
for k = 1:numel(Nvals)
    N = Nvals(k);
    
    % expRK4
    expRK4_approx(k) = expRK4(f, interval, u0, N);

    % sspRK3
    sspRK3_approx(k) = sspRK3(f, interval, u0, N);

    % sspRK2
    sspRK2_approx(k) = sspRK2(f, interval, u0, N);

end


% compute errors
expRK4_errors = abs(exact - expRK4_approx);
sspRK3_errors = abs(exact - sspRK3_approx);
sspRK2_errors = abs(exact - sspRK2_approx);

% display order table
expRK4_order = log2(expRK4_errors(1:end-1) ./ expRK4_errors(2:end));
sspRK3_order = log2(sspRK3_errors(1:end-1) ./ sspRK3_errors(2:end));
sspRK2_order = log2(sspRK2_errors(1:end-1) ./ sspRK2_errors(2:end));
header = {'expRK4', 'sspRK3', 'sspRK2'};
order_table = [expRK4_order, sspRK3_order, sspRK2_order];
disp('Order table');
disp(array2table(order_table, 'VariableNames',header));

% graph
deltaXVals = ((interval(2) - interval(1)) ./ Nvals)';

loglog(deltaXVals, expRK4_errors, 'b-', 'linewidth', 1.5); hold on;
loglog(deltaXVals, deltaXVals .^ (4), 'b-.', 'linewidth', 1.5);

loglog(deltaXVals, sspRK3_errors, 'r-', 'linewidth', 1.5);
loglog(deltaXVals, deltaXVals .^ (3), 'r-.', 'linewidth', 1.5);

loglog(deltaXVals, sspRK2_errors, 'g-', 'linewidth', 1.5);
loglog(deltaXVals, deltaXVals .^ (2), 'g-.', 'linewidth', 1.5);


xlabel('(b - a) / N'); ylabel('Error'); title('Runge-Kutta Approximations');
legend('expRK4', 'Order 4', 'sspRK3', 'Order 3', 'sspRK2', 'Order 2', 'Location', 'northwest');









