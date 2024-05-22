% Runge-Kutta method combined with centered-difference method
% u_t = u_xx, 0<x<1, t>0
% u(0,t) = u(1,0) = 0
% u(x,0) = sin(2*pi*x)
% Centered differences + f. (or b.) Euler


%% Crank-Nicholson
clc; clear variables; close all

% Temporal convergence

Nx = 100;
xvals = linspace(0, 1, Nx + 1)' ;
dx = xvals(2) - xvals(1);

lambdavals = 0.1:0.02:0.5;
lambdavals = 0.1:0.02:10;
errors = zeros(numel(lambdavals), 1);

tf = 0.3;
dtvals = lambdavals .* (dx);

Dxx = (1/(dx^2))*gallery('tridiag', Nx - 1, 1, -2, 1);    

for k = 1:numel(lambdavals)
    tvals = 0:dtvals(k):tf;
    if tvals(end) ~= tf
        tvals = [tvals, tf];
    end
    Nt = numel(tvals);

    u = sin(pi*xvals);
    u = u(2:end-1); % only on interior nodes (Dirichlet BCs)


    for i = 2:Nt
        dt = tvals(i) - tvals(i-1);
        u = (eye(Nx - 1) - (dt*Dxx/2)) \ ((eye(Nx - 1) + (dt*Dxx/2)) * u);
    end
    
    u = [0;u;0]; %bring BCs back in
    exact = exp(-pi^2*tf)*sin(pi*xvals);
    errors(k) = dx*sum(abs(u - exact)); %L1 norm
end

figure(1); clf; loglog(lambdavals, errors, 'b-', 'LineWidth', 1.5);
hold on; loglog(lambdavals, lambdavals.^2, 'b-.', 'LineWidth', 1.5);

% Spatial Convergence
Nxvals = [40, 80, 160, 320];
errors = zeros(numel(Nxvals), 1);
tf = 0.3;

for k = 1:numel(Nxvals)
    Nx = Nxvals(k);
    xvals = linspace(0, 1, Nx + 1)';

    dx = xvals(2) - xvals(1);
    Dxx = (1/(dx^2))*gallery('tridiag', Nx - 1, 1, -2, 1);  % centered difference method 
    
    dt = dx;    
    tvals = 0:dt:tf;

    if tvals(end) ~= tf
        tvals = [tvals, tf];
    end
    Nt = numel(tvals);

    u = sin(pi*xvals);
    u = u(2:end-1); % only on interior nodes (Dirichlet BCs)

    for i = 2:Nt
        dt = tvals(i) - tvals(i-1);
        u = (eye(Nx - 1) - (dt*Dxx/2)) \ ((eye(Nx - 1) + (dt*Dxx/2)) * u);
    end
    
    u = [0;u;0]; %bring BCs back in
    exact = exp(-pi^2*tf)*sin(pi*xvals);
    errors(k) = dx*sum(abs(u - exact)); % L1 norm
end

% Spatial accuracy
disp('Spatial Accuracy Order = ')
disp(log2(errors(1:end-1)./errors(2:end)))

%% RK 4
clc; clear variables;

Nx = 100;

% butcher table
A = [0, 0, 0, 0;
     0.5, 0, 0, 0;
     0, 0.5, 0, 0;
     0, 0, 1, 0];
b = [1/6, 1/3, 1/3, 1/6]';
c = [0, 0.5, 0.5, 1]';

s = numel(b); % length of butcher table


% Temporal convergence

Nx = 100;
xvals = linspace(0, 1, Nx + 1)' ;
dx = xvals(2) - xvals(1);

lambdavals = 0.1:0.02:0.5;
errors = zeros(numel(lambdavals), 1);

tf = 0.3;
dtvals = lambdavals .* (dx^2);

Dxx = (1/(dx^2))*gallery('tridiag', Nx - 1, 1, -2, 1);    

for k = 1:numel(lambdavals)
    tvals = 0:dtvals(k):tf;
    if tvals(end) ~= tf
        tvals = [tvals, tf];
    end
    Nt = numel(tvals);

    u = sin(pi*xvals);
    u = u(2:end-1); % only on interior nodes (Dirichlet BCs)

    % RK4 Method ---->
    % loop over timesteps
    kvals = zeros(s, Nx - 1); % initialize k vector
    for n = 2:Nt
        dt = tvals(n) - tvals(n-1);

        % for each partition, compute approximation of u
        for i = 1:s
            kvals(i, :) = Dxx*(u + (dt*(sum( A(i, :)' .* kvals)')));
        end
        % 
        % k1 = Dxx*u;
        % k2 = Dxx*(u+(dt*(0.5 .* k1)));
        % k3 = Dxx*(u+(dt*(0.5 .* k2)));
        % k4 = Dxx*(u+(dt*( k3)));
    
        % take weighted average of k values to approximate u
        u = u + dt*sum(b .* kvals)';
        % u = u + (dt/6)*(k1+(2*k2)+(2*k3)+k4);
    end
    u = [0;u;0]; %bring BCs back in
    exact = exp(-pi^2*tf)*sin(pi*xvals);
    errors(k) = dx*sum(abs(u - exact)); % L1 norm
end

figure(1); clf; loglog(lambdavals, errors, 'b-', 'LineWidth', 1.5);
hold on; loglog(lambdavals, lambdavals.^4, 'b-.', 'LineWidth', 1.5);

% Spatial Convergence
Nxvals = [40, 80, 160, 320];
errors = zeros(numel(Nxvals), 1);
tf = 0.3;

for k = 1:numel(Nxvals)
    Nx = Nxvals(k);
    xvals = linspace(0, 1, Nx + 1)';

    dx = xvals(2) - xvals(1);
    Dxx = (1/(dx^2))*gallery('tridiag', Nx - 1, 1, -2, 1);  % centered difference method 
    
    dt = 0.1*dx.^2;    
    tvals = 0:dt:tf;

    if tvals(end) ~= tf
        tvals = [tvals, tf];
    end
    Nt = numel(tvals);

    u = sin(pi*xvals);
    u = u(2:end-1); % only on interior nodes (Dirichlet BCs)

    % RK4 Method ---->
    % loop over timesteps
    kvals = zeros(s, Nx - 1); % initialize k vector
    for n = 2:Nt
        dt = tvals(n) - tvals(n-1);

        %for each partition, compute approximation of u
        for i = 1:s
            kvals(i, :) = (Dxx*u + (dt*(sum( A(i, :)' .* kvals)')));
        end
    
        
        % take weighted average of k values to approximate u
        u = u + dt*sum(b .* kvals)';
    end
    u = [0;u;0]; %bring BCs back in
    exact = exp(-pi^2*tf)*sin(pi*xvals);
    errors(k) = dx*sum(abs(u - exact)); % L1 norm
end

% Spatial accuracy
disp('Spatial Accuracy Order = ')
disp(log2(errors(1:end-1)./errors(2:end)))



















