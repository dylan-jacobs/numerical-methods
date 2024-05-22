clc; clear variables; close all;

Nx = 300;
xvals = linspace(0, 2*pi, Nx + 1)'; % Nx + 1 points
xvals = xvals(1:end-1); % chop off end point
dx = xvals(2) - xvals(1);
tf = 0.3;

Dxx = gallery('tridiag', Nx, 1, -2, 1);
Dx = gallery('tridiag', Nx, -1, 1, 0);

Dxx(1, end) = 1; Dxx(end, 1) = 1; % Periodic BCs
Dx(1, end) = -1;

Dxx = (1/(dx.^2))*Dxx;
Dx = (1/dx)*Dx;

lambdavals = (0.1:0.02:0.5)';
exact = exp(-tf).*sin(xvals-tf);
lambdaLen = numel(lambdavals);
errors = zeros(lambdaLen, 1);

for k = 1:lambdaLen
    dt = lambdavals(k).*dx.^2;
    tvals = (0:dt:tf)';

    if (tvals(end) ~= tf)
        tvals = [tvals; tf];
    end

    u = sin(xvals); % IC
    for n = 2:numel(tvals)
        dt = tvals(n) - tvals(n-1);
        u = (eye(Nx) - (dt*Dxx))\((eye(Nx) - (dt*Dx))*u);
    end

    errors(k) = dx*sum(abs(exact - u));

end

% Temporal error plot
figure(1);
loglog(lambdavals, errors, 'r-'); hold on;
loglog(lambdavals, lambdavals, 'r--');
title('Upwind error plot');
legend('Upwind', 'Order 1');


% SPATIAL ERROR
Nxvals = [40, 80, 160, 320];
errors = zeros(numel(Nxvals), 1);

for k = 1:numel(Nxvals)
    Nx = Nxvals(k);
    xvals = linspace(0, 2*pi, Nx + 1)';
    xvals = xvals(1:end-1);
    exact = exp(-tf).*sin(xvals-tf);
    
    dx = xvals(2) - xvals(1);
    
    Dxx = gallery('tridiag', Nx, 1, -2, 1);
    Dx = gallery('tridiag', Nx, -1, 1, 0);
    
    Dxx(1, end) = 1; Dxx(end, 1) = 1; % Periodic BCs
    Dx(1, end) = -1;
    
    Dxx = (1/(dx.^2))*Dxx;
    Dx = (1/dx)*Dx;

    dt = 0.1*dx.^2;
    tvals = 0:dt:tf;
    if (tvals(end) ~= tf)
        tvals(end) = tf;
    end

    u = sin(xvals); % IC
    for n = 2:numel(tvals)
        dt = tvals(n) - tvals(n-1);
        u = (eye(Nx) - (dt*Dxx))\((eye(Nx) - (dt*Dx))*u);
    end

    errors(k) = dx*sum(abs(exact - u));
end

% Spatial accuracy
disp('Upwind Spatial Accuracy Order = ')
disp(log2(errors(1:end-1)./errors(2:end)))















