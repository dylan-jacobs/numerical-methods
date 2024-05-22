% MATLAB Exercises
%% Meeting 1: 1/23/2024
clear variables; % clears vars
close all; % loses figures
clc; % clears command window

%% For loops

for i = 1:10
    disp( i*cos(i) ) % disp displays value; without disp, assigns var to 'ans'
end

%% while loops
n = 0;
term = 1;
tolerance = 0.05 * exp(2);
approx = 0;

while term > tolerance
    approx = approx + term;
    n = n + 1;
    term = 2^n / factorial(n);
end
disp('error')
disp(exp(2) - approx)
disp(['absolute error = ', num2str(exp(2) - approx)])


%% If statements
x = 5;
y = -1;

if x+y == 4
    disp(x*y);
end

if x > 0 && y > 0
    disp(x*y);
elseif x < 0 && y > 0
    disp(x+y);
else
    disp(0)
end

%% Vectors and matrices

% % % HOW TO BUILD VECTORS AND MATRICES
% Vectors:
vec1 = 1:2:7; % default vector orientation is row-vector
vec1 = vec1'; % single quote after variable takes transpose
vec2 = linspace(0,2,5)'; % linspace(start, end, num of points)
mat1 = [1,2,3 ; 4,5,6 ; 7,8,9]; % construct vector/matrix using brackets (commas separate rows, semicolons separate rows)

% Example: 5x1 vector with entries unknown
    % Method 1: appending existing array --> VERY INEFFICIENT!
    vec3 = [];
    for k = 1:5
        vec3 = [vec3; 3*k];
    end

    % Method 2: preallocating array size
    vec4 = zeros(5,1);
    for k = 1:5
        vec4(k) = 3*k;
    end

% Example: 5x5 matrix; (i,j) entry is i+j
    % Method 1: nested for loops --> VERY INEFFICIENT!
    mat2 = zeros(5,5);
    for i = 1:5
        for j = 1:5
            mat2(i,j) = i+j;
        end
    end

    % Method 2: Vectorize it --> More efficient 
    mat3 = (1:5) + (1:5)'; % column vector + row vector = added matrix
    % Equivalent command: repmat(array, # rows, # cols)
    
%% Matrix, vector multiplication
clear variables;
close all;

A = [1,2,3 ; 4,5,6 ; 7,8,9];
B = [1,0,0 ; 1,1,1 ; 0,0,1];
x = [ 1,1,2 ]';

% Standard matrix multiplication:
standard_product = A * B;

% Component-wise multiplication:
component_product = A .* B;

% Example: Matrix-vector multiplication
standard_mtx_vector_product = A*x;
component_wise_mtx_vector_product = A.*x; % --> MATLAB autofills matrix if dimensions don't match between vectors, matrices























