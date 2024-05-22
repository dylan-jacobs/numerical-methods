%% Numerical integration methods
 clc;
 clear variables;
 close all;


 %% Single run
 a = 0; b = 2; % [a, b]
 N = 10; % N intervals, N+1 points
 f = @(x) sin(x);

exact = (-cos(b)) + cos(a);
approx = left_reimman(a, b, N, f);
error = abs(exact - approx);

disp(error);

%% Testing Order of Accuracy
% Determine error from graph
% E ~ deltaX^P
% log(E) ~ P * log(deltaX)
% P = log(E) / log(deltaX) 
% P = slope of log(E) plotted over log(deltaX)

a = 0; b = 2; % [a, b]
f = @(x) sin(x);
Nvals= [4, 8, 16, 32, 64, 128, 256]; % N intervals, N+1 points

exact = (-cos(b)) + cos(a);
errors = zeros(numel(Nvals), 1);

for k = 1:numel(Nvals)
    N = Nvals(k);
    approx = right_reimman(a, b, N, f);
    error = abs(exact - approx);
    errors(k) = error;
end

% Plot the error
figure(1); clf;
loglog((b-a)./Nvals, errors, 'b-', 'linewidth', 1.5);
xlabel('(b-a)/N'); ylabel('Absolute error'); title('Composite Midpoint Rule');

figure(1);hold on;loglog((b-a)./Nvals, ((b-a)./Nvals).^1, 'r-', 'linewidth', 1.5);
legend('midpoint', 'order 2', 'location', 'northwest');

%% HOMEWORK
% Left Reimman Sum, Right Reimman Sum (1st order), Trapezoid Rule (2nd order)

% Maybe rewrite order of accuracy test script if time
% Determine error from graph
% E ~ deltaX^P
% log(E) ~ P * log(deltaX)
% P = log(E) / log(deltaX) 
% P = slope of log(E) plotted over log(deltaX)

a = 0; b = 2; % interval [a, b]
Nvals = [4, 8, 16, 32, 64, 128, 256]; % N intervals
f = @(x) sin(x);
exact_integral = -cos(b) + cos(a);

errors = zeros(4, numel(Nvals));

for i = 1:numel(Nvals)
    N = Nvals(i);
    errors(1, i) = left_reimman(a, b, N, f);
    errors(2, i) = right_reimman(a, b, N, f);
    errors(3, i) = midpoint(a, b, N, f);
    errors(4, i) = trapezoid_sum(a, b, N, f);
end

errors = abs(errors - exact_integral);

deltaXVals = (b-a) ./ Nvals;
figure(1); clf; loglog(deltaXVals, errors(1, :), 'b-', 'linewidth', 1.5);
figure(1); hold on; loglog(deltaXVals, errors(2, :), 'r-', 'linewidth', 1.5);
figure(1); hold on; loglog(deltaXVals, errors(3, :), 'g-', 'linewidth', 1.5);
figure(1); hold on; loglog(deltaXVals, errors(4, :), 'm-', 'linewidth', 1.5);

% Graph Order 2 accuracy function
figure(1); hold on; loglog(deltaXVals, deltaXVals .^ 2, 'r--', 'linewidth', 1.5);

% Graph Order 1 accuracy function
figure(1); hold on; loglog(deltaXVals, deltaXVals, 'b--', 'linewidth', 1.5);

legend('Left Reimman approximation', 'Right Reimman approximation', 'Midpoint approximation', 'Trapezoid sum approximation', 'Order 2 function', 'Order 1 function', 'location', 'southeast');
xlabel('(a-b)/N'); ylabel('Error'); title('Order of Accuracy as a function of deltaX')










































 