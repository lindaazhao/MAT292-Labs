function [IEM_t, IEM_y] = IEM_solver(f, t0, tN, y0, h)
% f = function (right side of ODE)
% [t0, tN] = internal on which to solve the ODE
% y0 = initial condition of the ODE
% h = step size
    if mod(tN-t0, h) ~= 0                           % Range does not fit an integer multiple of the step size
        tN = t0 + (floor((tN-t0)/h) + 1) * h;       % Increase tN such that it now fits an integer multiple of h
    end
    
    IEM_t = t0:h:tN;
    
    IEM_y = NaN(1, length(IEM_t));
    IEM_y(1) = y0;

    for j = 2:length(IEM_t)-1
        IEM_y(j) = IEM_y(j-1) + h*f(IEM_t(j-1), IEM_y(j-1));
        IEM_y(j+1) = IEM_y(j) + h*f(IEM_t(j), IEM_y(j));
        IEM_y(j) = IEM_y(j-1) + (h/2)*(f(IEM_t(j-1), IEM_y(j-1)) + f(IEM_t(j), IEM_y(j)));
    end
end