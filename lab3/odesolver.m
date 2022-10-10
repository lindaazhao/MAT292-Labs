function [odesoln] = odesolver(f1, t0, tN, y0, h)
% f1 = function (right side of ODE)
% [t0, tN] = internal on which to solve the ODE
% y0 = initial condition of the ODE
% h = step size
    if mod(tN-t0, h) ~= 0                           % Range does not fit an integer multiple of the step size
        tN = t0 + (floor((tN-t0)/h) + 1) * h;       % Increase tN such that it now fits an integer multiple of h
    end
    
    soln_t = t0:h:tN;
    
    odesoln = NaN(1, length(soln_t));
    odesoln(1) = y0;

    for j = 2:length(soln_t)-1
        odesoln(j) = odesoln(j-1) + h*f1(soln_t(j-1), odesoln(j-1));
        odesoln(j+1) = odesoln(j) + h*f1(soln_t(j), odesoln(j));
        odesoln(j) = odesoln(j-1) + (h/2)*(f1(soln_t(j-1), odesoln(j-1)) + f1(soln_t(j), odesoln(j)));
    end
end
