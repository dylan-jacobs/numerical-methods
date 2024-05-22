% approximate function f using sspRK3 method
% u0 = initial condition
% N points

function [approx] = sspRK3(f, interval, u0, N)

    % butcher table
    A = [0, 0, 0;
         1, 0, 0;
         0.25, 0.25, 0];
    b = [1/6, 1/6, 2/3]';
    c = [0, 1, 0.5]';
    
    approx = runge_kutta(A, b, c, f, u0, interval, N);

end