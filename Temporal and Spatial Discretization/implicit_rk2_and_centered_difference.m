%% RK 2
clc; clear variables;

% butcher table
gamma = 1 - (1/sqrt(2));
A = [gamma, 0;
     1-(2*gamma), gamma];
b = [1/2, 1/2]';
c = [gamma, 1-gamma]';

s = numel(b); % length of butcher table


% Temporal convergence
Nx = 100;
xvals = linspace(0, 1, Nx + 1)' ;
dx = xvals(2) - xvals(1);

lambdavals = 0.05:0.05:10;
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

    % RK4 Method ---->
    % loop over timesteps
    kvals = zeros(s, Nx - 1); % initialize k vector
    for n = 2:Nt
        dt = tvals(n) - tvals(n-1);

        k1 = (eye(Nx - 1) - (Dxx*dt*A(1, 1)))\(Dxx*u);
        k2 = (eye(Nx - 1) - (Dxx*dt*A(2, 2)))\(Dxx*u + Dxx*dt*(A(2, 1)*k1));
    
        % take weighted average of k values to approximate u
        u = u + dt*sum(b .* [k1, k2]')';
    end
    u = [0;u;0]; %bring BCs back in
    exact = exp(-pi^2*tf)*sin(pi*xvals);
    errors(k) = dx*sum(abs(u - exact)); % L1 norm
end

figure(1); clf; loglog(lambdavals, errors, 'b-', 'LineWidth', 1.5);
hold on; loglog(lambdavals, lambdavals.^2, 'b-.', 'LineWidth', 1.5);
xlabel('\lambda'); ylabel('Error');
legend('RK2', 'Order 2');

% Spatial Convergence
Nxvals = [40, 80, 160, 320];
errors = zeros(numel(Nxvals), 1);
tf = 0.3;

for k = 1:numel(Nxvals)
    Nx = Nxvals(k);
    xvals = linspace(0, 1, Nx + 1)';

    dx = xvals(2) - xvals(1);
    Dxx = (1/(dx^2))*gallery('tridiag', Nx - 1, 1, -2, 1);  % centered difference method 
    
    dt = 1*dx;    
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
        
        k1 = (eye(Nx - 1) - (Dxx*dt*A(1, 1)))\(Dxx*u);
        k2 = (eye(Nx - 1) - (Dxx*dt*A(2, 2)))\(Dxx*u + Dxx*dt*(A(2, 1)*k1));
    
        % take weighted average of k values to approximate u
        u = u + dt*sum(b .* [k1, k2]')';
            
    end
    u = [0;u;0]; %bring BCs back in
    exact = exp(-pi^2*tf)*sin(pi*xvals);
    errors(k) = dx*sum(abs(u - exact)); % L1 norm
end

% Spatial accuracy
disp('Spatial Accuracy Order = ')
disp(log2(errors(1:end-1)./errors(2:end)))