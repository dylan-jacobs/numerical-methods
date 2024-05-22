
%% DOWNWIND
clc; clear variables; close all;

tf = 1;
u0 = @(x) sin(x);

Nx = 100;
xvals = linspace(0, 2*pi, Nx + 1)';
dx = xvals(2) - xvals(1);

a = 1;
Dx = gallery('tridiag', Nx+1, 0, -1, 1);
Dx(end, 1) = 1;

lambdavals = (0.1:0.02:0.5)';

errors = zeros(numel(lambdavals), 1);
exact = sin(xvals - tf);
 
for k = 1:numel(lambdavals)
    dt = lambdavals(k)*dx;
    tvals = (0:dt:tf)';

    u = u0(xvals);

    if (tvals(end) ~= tf)
        tvals = [tvals; tf];
    end

    for i = 2:numel(tvals)
        u = u - (a*dt/dx).*(Dx*u);
    end

    u = u(1:end-1);
    errors(k) = dx*sum(abs(exact(1:end-1) - u));
end

% Solution curves
figure(1);
plot(xvals(1:end-1), u); hold on;
plot(xvals(1:end-1), exact(1:end-1));
title('Downwind solution curves');
xlabel('x'); ylabel('u');
legend('approx', 'exact');

% Spatial Error
clear variables;
Nxvals = [40, 80, 160, 320];
a = 1;
u0 = @(x) sin(x);
tf = 1;

errors = zeros(numel(Nxvals), 1);

for k = 1:numel(Nxvals)
    Nx = Nxvals(k);
    xvals = linspace(0, 2*pi, Nx + 1)';
    dx = xvals(2) - xvals(1);
    dt = 0.1*dx;
    tvals = 0:dt:tf;
    if (tvals(end) ~= tf)
        tvals(end) = tf;
    end

    u = u0(xvals);
    Dx = gallery('tridiag', Nx+1, 0, -1, 1);
    Dx(end, 1) = 1;
    
    for i = 2:numel(tvals)
        dt = tvals(i)-tvals(i-1);
        u = u - (a*dt/dx).*(Dx*u);
    end
    u = u(1:end-1);
    exact = sin(xvals - tf);
    errors(k) = dx*sum(abs(u - exact(1:end-1))); % L1 norm
end

% Spatial accuracy
disp('Downwind Spatial Accuracy Order = ')
disp(log2(errors(1:end-1)./errors(2:end)))







%% UPWIND

clc; clear variables; close all;

tf = 1;
u0 = @(x) sin(x);

Nx = 800;
xvals = linspace(0, 2*pi, Nx + 1)';
dx = xvals(2) - xvals(1);

a = 1;
Dx = gallery('tridiag', Nx+1, -1, 1, 0);
Dx(1, end) = -1;

lambdavals = (0.1:0.02:0.5)';

errors = zeros(numel(lambdavals), 1);
exact = sin(xvals - tf);
 
for k = 1:numel(lambdavals)
    dt = lambdavals(k)*dx;
    tvals = (0:dt:tf)';

    u = u0(xvals);

    if (tvals(end) ~= tf)
        tvals = [tvals; tf];
    end

    for i = 2:numel(tvals)
        u = u - (a*dt/dx).*(Dx*u);
    end

    u = u(1:end-1);
    errors(k) = dx*sum(abs(exact(1:end-1) - u));
end

% Temporal error plot
figure(1);
loglog(lambdavals, errors, 'r-'); hold on;
loglog(lambdavals, lambdavals, 'r--');
title('Upwind error plot');
legend('Upwind', 'Order 1');

% Solution curves
figure(2);
plot(xvals(1:end-1), u); hold on;
plot(xvals(1:end-1), exact(1:end-1));
title('Upwind solution curves');
xlabel('x'); ylabel('u');
legend('approx', 'exact');

% Spatial Error
clear variables;
Nxvals = [40, 80, 160, 320];
a = 1;
u0 = @(x) sin(x);
tf = 1;

errors = zeros(numel(Nxvals), 1);

for k = 1:numel(Nxvals)
    Nx = Nxvals(k);
    xvals = linspace(0, 2*pi, Nx + 1)';
    dx = xvals(2) - xvals(1);
    dt = 0.1*dx;
    tvals = 0:dt:tf;
    if (tvals(end) ~= tf)
        tvals(end) = tf;
    end

    u = u0(xvals);
    Dx = gallery('tridiag', Nx+1, -1, 1, 0);
    Dx(1, end) = -1;
    
    for i = 2:numel(tvals)
        dt = tvals(i)-tvals(i-1);
        u = u - (a*dt/dx).*(Dx*u);
    end
    u = u(1:end-1);
    exact = sin(xvals - tf);
    errors(k) = dx*sum(abs(u - exact(1:end-1))); % L1 norm
end

% Spatial accuracy
disp('Upwind Spatial Accuracy Order = ')
disp(log2(errors(1:end-1)./errors(2:end)))











