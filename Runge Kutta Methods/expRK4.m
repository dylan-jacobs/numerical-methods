% approximate function f using expRK4 method
% u0 = initial condition
% N points

function [approx] = expRK4(f, interval, u0, N)

    % butcher table
    A = [0, 0, 0, 0;
         0.5, 0, 0, 0;
         0, 0.5, 0, 0;
         0, 0, 1, 0];
    b = [1/6, 1/3, 1/3, 1/6]';
    c = [0, 0.5, 0.5, 1]';
    
    approx = runge_kutta(A, b, c, f, u0, interval, N);

end