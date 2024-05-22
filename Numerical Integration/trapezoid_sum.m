% left Reimman sum approximation

function output = trapezoid_sum(a, b, N, f)
    % interval [a, b]
    % N intervals
    % function f
    
    xvals = linspace(a, b, N+1);
    deltaX = xvals(2) - xvals(1);
    output = 0;

    for i = 1:(numel(xvals) - 1)
        output = output + ((f(xvals(i)) + f(xvals(i+1))) / 2);
    end
    output = output * deltaX;
end