% approximate function f using sspRK2 method
% u0 = initial condition
% N points

function [approx] = sspRK2(f, interval, u0, N)

    gamma = 1-(sqrt(2)/2);
    delta = 1-(1/(2*gamma));

    % butcher table
    A = [0, 0, 0;
         gamma, 0, 0;
         delta, 1-delta, 0];
    b = [delta, 1-delta, 0]';
    c = [0, gamma, 1]';
    
    approx = runge_kutta(A, b, c, f, u0, interval, N);

end