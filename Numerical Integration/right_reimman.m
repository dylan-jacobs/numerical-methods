% right Reimman sum approximation

function output = right_reimman(a, b, N, f)
    % interval [a, b]
    % N intervals
    % function f
    
    xvals = linspace(a, b, N+1);
    deltaX = xvals(2) - xvals(1);
    output = 0;

    for i = 2:numel(xvals)
        output = output + f(xvals(i));
    end
    output = output * deltaX;
end