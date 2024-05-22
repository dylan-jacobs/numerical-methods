clear all; close all;

%% Accuracy in time

Nx = 100;
xvals = linspace(0,1,Nx+1)'; %Nx+1 points
dx = xvals(2)-xvals(1);

%lambdavals = 0.1:0.02:10; %dt = lambda*dx for b. Euler
lambdavals = 0.01:0.02:0.5; %dt = lambda*dx^2 for f. Euler
errvals = zeros(numel(lambdavals),1);

Tf = 0.3;

for k = 1:numel(lambdavals)
    dt = lambdavals(k)*dx^2
    tvals = 0:dt:Tf;
    if tvals(end) ~= Tf
        tvals = [tvals, Tf];
    end
    Nt = numel(tvals);

    u = sin(pi*xvals); %IC

    u = u(2:end-1); %only on interior nodes bc Dirichlet BCs
    
    Dxx = (1/dx^2)*gallery('tridiag',Nx-1,1,-2,1);

    for n = 2:Nt
        dtn = tvals(n) - tvals(n-1);
        %u = (eye(Nx-1) - dtn*Dxx)\u; %b. Euler
        u = u + dtn*Dxx*u; %f. Euler
    end
    u = [0;u;0]; %bring BCs back in
    exact = exp(-pi^2*Tf)*sin(pi*xvals);
    errvals(k) = dx*sum(abs(u - exact)); %L1 norm
end


figure(1);clf;loglog(lambdavals,errvals,'b-','linewidth',1.5);
hold on;loglog(lambdavals,0.001*lambdavals.^1,'black-.','linewidth',1.5);
xlabel('\lambda');ylabel('L1 error');


%% Convergence under spatial mesh refinement
clear all;

Nxvals = [40,80,160,320];
errvals = zeros(numel(Nxvals),1);

for k = 1:numel(Nxvals)
    Nx = Nxvals(k);
    xvals = linspace(0,1,Nx+1)'; %Nx+1 points
    dx = xvals(2)-xvals(1);
    
    lambda = 3;%0.3; %dt = lambda*dx^2 for f. Euler
    dt = lambda*dx;
    Tf = 0.3;
    tvals = 0:dt:Tf;
    if tvals(end) ~= Tf
        tvals = [tvals, Tf];
    end
    Nt = numel(tvals);

    u = sin(pi*xvals); %IC

    u = u(2:end-1); %only on interior nodes bc Dirichlet BCs
    
    Dxx = (1/dx^2)*gallery('tridiag',Nx-1,1,-2,1);

    for n = 2:Nt
        dtn = tvals(n) - tvals(n-1);
        u = (eye(Nx-1) - dtn*Dxx)\u; %b. Euler
        %u = u + dtn*Dxx*u; %f. Euler
    end
    u = [0;u;0]; %bring BCs back in
    exact = exp(-pi^2*Tf)*sin(pi*xvals);
    errvals(k) = dx*sum(abs(u - exact)); %L1 norm
end

disp('Order = ')
disp(log2(errvals(1:end-1)./errvals(2:end)))



