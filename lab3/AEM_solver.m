function [AEM_t, AEM_y] = AEM_solver(f, t0, tN, y0, h)
% f = function (right side of ODE)
% [t0, tN] = internal on which to solve the ODE
% y0 = initial condition of the ODE
% h = initial step size
    
    AEM_t = [t0];
    AEM_y = [y0];
    tol = 1e-8;
    j = 2; % Start the loop by calculating the second value for AEM_y

    while (AEM_t(end)<tN)
        while true
            Y = AEM_y(j-1) + h*f(AEM_t(j-1), AEM_y(j-1));               % Take one Euler step size of h
            Z_half = AEM_y(j-1) + (h/2)*f(AEM_t(j-1), AEM_y(j-1));      % Take first half Euler step size of h/2
            Z = Z_half + (h/2)*f(AEM_t(j-1)+(h/2), Z_half);             % Take second half Euler step size of h/2

            D = Z-Y;

            if abs(D) < tol                         % Step was successful
                AEM_y = [AEM_y (Z + D)];            % AEM_y will have local error of O(h^3)
                AEM_t = [AEM_t (AEM_t(j-1)+h)];
                j = j + 1;
                break
            else
                h = 0.9*h*min(max(tol/abs(D),0.3),2);
            end
        end
    end
end

