function [T,X] = solvesystem_zhaoli50(f, g, t0, tN, x0, h)
% f,g = right side of ODE functions
% [t0, tN] = interval on which to solve the ODE
% x0 = vector containing initial condition of the ODE (x0, y0)
% h = step size
    if mod(tN-t0, h) ~= 0                               % Range does not fit an integer multiple of the step size
            tN = t0 + (floor((tN-t0)/h) + 1) * h;       % Increase tN such that it now fits an integer multiple of h
    end
        
    T = t0:h:tN;
    
    x1 = zeros(1, length(T));
    x2 = zeros(1, length(T));
    
    x1(1) = x0(1);
    x2(1) = x0(2);

    for j = 2:length(T)-1
        x1(j) = x1(j-1) + h*f(T(j-1), x1(j-1), x2(j-1));                                                   % Euler step
        x1(j) = x1(j-1) + (h/2)*(f(T(j-1), x1(j-1), x2(j-1)) + f(T(j), x1(j), x2(j)));    % Improved Euler step

        x2(j) = x2(j-1) + h*g(T(j-1), x1(j-1), x2(j-1));                                                   % Euler step
        x2(j) = x2(j-1) + (h/2)*(g(T(j-1), x1(j-1), x2(j-1)) + g(T(j), x1(j), x2(j)));
    end

    X = NaN(2, length(T));
    X(1,:) = x1;
    X(2,:) = x2;
end

