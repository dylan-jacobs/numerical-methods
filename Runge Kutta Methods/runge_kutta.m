% Runge-Kutta Method
% A, b, c compose the butcher table
% function f
% u0 = initial condition
% N points

function [u] = runge_kutta(A, b, c, f, u0, interval, N)

    s = numel(b); % length of butcher table
    tvals = linspace(interval(1), interval(2), N); % partition interval into N timesteps
    dt = (tvals(2)-tvals(1)); % delta t
    kvals = zeros(s, 1); % initialize k vector

    % initialize u
    u = u0;

    % loop over partitions
    for n = 2:N

        % for each partition, compute approximation of u
        for i = 1:s
            kvals(i) = f(u + dt*(sum( A(i, :)' .* kvals)), tvals(n-1) + dt*c(i));
        end
        
        % take weighted average of k values to approximate u
        u = u + dt*sum(b .* kvals);

    end

end