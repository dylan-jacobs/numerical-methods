clear variables;
clc;

%% If Statements: Problem 2
x = -5;
y = 1;

if x>0
    if y>0
        z = x + 2;
    else
        z = x - 2;
    end
else
    if y < 0
        z = y + 2;
    else 
        z = y - 2;
    end 
end

% Z = -1

%% Loops: Problem 2
clear variables;
clc;

x = 0;
tol = 50.2;
while x < tol
    x = x + 1;
end
disp(['Final value of x: ', num2str(x)]);

% Final value of x: 51

%% Loops: Problem 3
clear variables;
clc;

y = 2;
x = y;
for k = 1:5
    y = x + k;
    x = y - 1;
end
disp(['Final value of x is ', num2str(x)]);
disp(['Final value of y is ', num2str(y)]);

% Final value of x is 12
% Final value of y is 13

%% Loops: Problem 4
clear variables;
clc;

x = input('Enter the value of x for the maclaurin approximation of e^x where -1 < x < 1.');
while ~ ((-1 < x) && (x < 1))
    x = input('Input the value of x for the maclaurin approximation of e^x where -1 < x < 1.');
end

n = input('Enter the number of terms to sum to approximate e^x: ');
while n < 1
    n = input('Enter the number of terms to sum to approximate e^x: ');
end

approx = 1;

for i = 1:n
    approx = approx + ((x^i) / factorial(i));
end

disp([sprintf('Actual value of e^%f ', x), num2str(exp(x))])
disp([sprintf('Approximation of e^%f = ', x), num2str(approx)]);


%% Loops: Problem 5
clear variables;
clc;

count = 0;
for i = 1:6
    for j = 1:6
        if i+j > 10
            count = count + 1;
        end
    end
end

fprintf('There are %d possible combinations of 2 6-sided die that are greater than 10', count);
fprintf('\n');

%% Matrices and Vectors: Problem 1
clear variables;
xvals = linspace(0, 2 * pi, 5);
yvals = linspace(0, pi, 3);

% For-loop method:
fxnvals1 = zeros(5,3);
for i = 1:5
    for j = 1:3
        fxnvals1(i,j) = sin(xvals(i) + yvals(j));
    end 
end
disp('Matrix resulting from for loop:');
disp(fxnvals1);

% Vector multiplication method
fxnvals2 = sin(xvals' + yvals);
disp('Matrix resulting from vectorization:');
disp(fxnvals2);

%% Matries and Vectors: Problem 3
clear variables;
clc;

xvals = linspace(0,1,11);
yvals = linspace(1,2,11);

% Consider the function f(x, y) = xy − 0.5y. Write a MATLAB®code that calculates the number
% of positive values when evaluating f(x, y) over all points (xvals(i),yvals(j)). There are many
% ways to do this.
% Answer: 55

solutions = (xvals' * yvals) - (0.5 * yvals);
numPositiveVals = solutions(solutions > 0);
disp(['The number of positive values in the solutions matrix of f(x) is ', num2str(size(numPositiveVals,1))]);


















