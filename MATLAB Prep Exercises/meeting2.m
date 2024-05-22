clear variables
close all

%% Functions Intro
% Example 1 (informal)
fxn_informal = @(x) (x^(1/3) + 1) * sin(x^2 + x + 1);

% Example 2
% We want to call f on multiple xvals ==> do them all at once in vector
f = @(x) (x.^(1/3) + 1) .* sin(x.^2 + x + 1);
xvals = (0:pi/4:pi)';
f(xvals);


%% Plotting
close all
xvals = linspace(-1*pi, pi, 100)';
f = @(x) sin(x);
g = @(x) cos(x);
yvals = f(xvals);
y2vals = g(xvals);

% Figure 1
figure(1);clf;plot(xvals, yvals, 'red', 'linewidth', 1.5);
xlabel('X');ylabel('y');title('Plot of sin(X)');
axis([-pi, pi, -1.2, 1.2]); % trims the graph ([xMin, xMax, yMin, yMax])

figure(2);clf;plot(xvals, yvals, 'red', 'linewidth', 1.5);
xlabel('X');ylabel('y');title('Plot of sin(X)');
axis([-pi, pi, -1.2, 1.2]); % trims the graph ([xMin, xMax, yMin, yMax])

% Figure 2
figure(2);hold on; plot(xvals, y2vals, 'blue', 'linewidth', 1.5);
legend('sin(x)', 'cos(x)', 'location', 'northeast')











%% Functions

% x = input, output = function value
function output = fxn1(x)
    output = (x^(1/3) + 1) * sin(x^2 + x + 1);
end

% x = input
% f = function value f(x)
% g = function value g(x)
function [f, g] = fxn2(x)
    f = sin(x);
    g = cos(x);
end






