function [T,Y] = DE2_zhaoli50(f, t0, tN, y0, y1, h)
% f = function y'' = f(t, y, y')
% [t0, tN] = interval on which to solve the ODE
% y0, y1 = initial conditions (y(t0) and y'(t0))
% h = step size

    if mod(tN-t0, h) ~= 0                               % Range does not fit an integer multiple of the step size
            tN = t0 + (floor((tN-t0)/h) + 1) * h;       % Increase tN such that it now fits an integer multiple of h
    end
        
    % Initialize t vector
    T = t0:h:tN;

    % Initialize solution vector
    Y = zeros(1, length(T));

    % Use initial conditions to get first two points
    Y(1) = y0;
    Y(2) = y0 + y1*h;

    for i = 2:length(T)-1
        % Approximate y' = Y_pr at i
        Y_pr = (Y(i) - Y(i-1))/h;

        % Evaluate y'' using t, y, and y' at i
        Y_dpr = f(T(i+1), Y(i), Y_pr);

        % Take next step for y
        Y(i+1) = (h^2)*Y_dpr + 2*Y(i)-Y(i-1);
    end
end

