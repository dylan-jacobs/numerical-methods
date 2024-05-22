clear variables;
clc;
close all;

%% Function Exercises 

% Exercise 1
% part a)
f = @(x) x^2 + x + 1;
disp(['f(2) = ', num2str(f(2))]);

% part b)
xvals = (0:0.2:1)';
for i = 1:numel(xvals)
    disp(['f(', num2str(xvals(i)), ') = ', num2str(f(xvals(i)))]);
end
% part c)
f2 = @(x) x .^ 2 + x + 1;
disp(mat2str(f2(xvals)));

% Exercise 2
a = [-1, 2, -1]';
b = [3, -5, 2]';
proj = dot(a, b) / norm(a);
disp(['The scalar projection of b onto a is ', num2str(proj)]);

% Exercise 3
f = @(x) (x .^ 2 + x - 1) ./ (x .^ 2 + 1);
xvec = linspace(0,1,10)';
disp(f(xvec));

%% Plotting Exercises
clear variables;
clc;
close all;

% Exercise 1
xvec = linspace(0, 20, 1000);
yvec = besselj(0,xvec);

figure(1); clf; plot(xvec, yvec, 'black-.', 'linewidth', 1.5);
xlabel('X'); ylabel('y'); title('Plotting the zeroth-order Bessel function');

% Exercise 2
X = linspace(-1, 1, 1000);
y1 = X.^2;
y2 = X.^3;
y3 = X.^4;

figure(2); hold on;
plot(X, X, 'b', 'linewidth', 1.5);
plot(X, y1, 'm', 'linewidth', 1.5);
plot(X, y2, 'r', 'linewidth', 1.5);
plot(X, y3, 'g', 'linewidth', 1.5);
xlabel('X');
ylabel('y');
title('Plotting various polynomials');
legend('X', 'X^2', 'X^3', 'X^4', Location='north');



















