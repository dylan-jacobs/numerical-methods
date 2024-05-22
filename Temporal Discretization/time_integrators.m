% Hard-coded example using forward/backward Euler Methods
% Goal: solve du/dt = f(u, t), where f(u, t) = -2tu, u(0) = 1
% u(t) = e^(-t^2))



clc; close all; clear variables;

Nvals = [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048];
errors = zeros(numel(Nvals), 1);
dtvals = zeros(numel(Nvals), 1); 
Tf = 2; % final time

for k = 1:numel(Nvals)
    N = Nvals(k);
    tvals = linspace(0, Tf, N + 1);
    dt = tvals(2) - tvals(1);
    dtvals(k) = dt;

    % initial condition
    u = 1; 
    for n = 2:N
        % Forward Euler Method
        % u = u - (dt*2*tvals(n - 1)*u);

        % Backward Euler Method
        % u = u / (1 + (2*dt*tvals(n)));

        % Crank-Nicholson Method
        u = u + ((dt / 2) * ((-2*tvals(n-1)*u) + (-2*tvals(n)*u)));
    end
    
    errors(k) = abs(exp(-Tf^2) - u);
end

log2(errors(1:end-1) ./ errors(2:end))
figure(1); clf; loglog(dtvals, errors, 'b-', 'LineWidth', 1.5);
hold on; loglog(dtvals, dtvals, 'b-.', 'LineWidth', 1.5);

































