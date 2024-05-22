% Composite midpoint rule 

function [output] = midpoint(a, b, N, f)
    % interval [a, b]
    % N intervals, N+1 points
    % function f
    xvals = linspace(a, b, N+1);
    dx = xvals(2) - xvals(1);
    
    mid_rule = 0; % midpoint rule approximation
    
    for i = 1:N
     mid_rule = mid_rule + dx * f((xvals(i) + xvals(i+1)) / 2);
    end
    
    output = mid_rule;
end
