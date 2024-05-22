% Create an error plot to determine a function's order

function [order] = CreateErrorPlot(a, b, exact_vals, approximations, NVals, expected_order, functionTitle)
errors = abs(approximations - exact_vals);
deltaXVals = (b - a) ./ NVals;

order = log2(errors(1:end-1) ./ errors(2:end));

figure(); clf; 
loglog(deltaXVals, errors, 'b-', 'linewidth', 1.5); hold on;
loglog(deltaXVals, deltaXVals .^ (expected_order), 'r-', 'linewidth', 1.5);
xlabel('(b - a) / N'); ylabel('Error'); title(sprintf('Composite Gauss-Legendre Quadrature Approximation for Function %s', functionTitle));
legend(sprintf('Gauss-Legendre Quadrature Approximation Error (M=%d)', expected_order), sprintf('Order-%d Approximation Function', expected_order), 'location', 'north');
hold off;