function [  ] = comp( a, n, l, dt, p0, w0, tol, dw )
%COMP Combination of phase map and Poincare map.
%   a - angle of attack (rad)
%   n - resolution of the phase map
%   l - number of solutions to compute
%   dt - time step of the solutions
%   p0 - starting roll angle (rad)
%   w0 - starting angular velocity (rad/tick)
%   tol - omega tolerance
%   dw - resolution of the p-map
    % Initialize stability derivatives
    tau = 0.5;
    c1 = 0.2;
    c2 = -0.1;
    c3 = 0.1;
    c4  = -1;
    c5 = -0.6;
    c6 = -0.02;
    % Calculate function values
    f = c1*sin(a)+c2*cos(a);
    g = sin(a)*(c3*(sin(a))^2+(c4+c2)*sin(a)*cos(a)-c1/2*(cos(a))^2);
    h = c6*(sin(a))^2-c5/6*(cos(a))^2;
    % Generate value matrices
    fval = linspace(-pi, pi, n);
    wval = linspace(-1, 1, n);
    [p, w] = meshgrid(fval, wval);
    dotp = NaN(size(p));
    dotw = NaN(size(w));
    % Fill the phase map
    for i = 1:numel(p)
        dotp(i) = w(i);
        dotw(i) = tau*(f+g*p(i)^2)*w(i)+sin(a)*(c5+h*p(i)^2)*p(i);
    end
    close all;
    figure;
    quiver(p, w, dotp, dotw, 'red');
    axis equal manual;
    xlabel('\phi');
    ylabel('\omega');
    title(['Phase map for \alpha = ' num2str(a)]);
    hold on;
    % Calculate trajectories using Euler's
    X = [[p0; w0] zeros(2,l-1)];
    for k = 1:l-1
        if((X(1,k) > pi) || (X(1,k) < -pi))
            X(1,k) = - X(1,k);
        end
        X(:,k+1) = X(:,k)+dt*[X(2,k);tau*(f+g*X(1,k)^2)*X(2,k)+sin(a)*(c5+h*X(1,k)^2)*X(1,k);];
    end
    scatter(p0,w0,'ko');
    scatter(X(1,l),X(2,l),'ks');
    scatter(X(1,:),X(2,:),'b.');
    hold off;
    % P-Map
    R = [];
    for w0 = 0:dw:1
        k = 1;
        X = [0; w0];
        while (~((X(1,k) < 0+tol) && (X(1,k) > 0-tol) && (X(2,k) > 0-tol)) || (k < (1/dt)))
            if((X(1,k) > pi) || (X(1,k) < -pi))
                X(1,k) = - X(1,k);
            end
            X(:,k+1) = X(:,k)+dt*[X(2,k);tau*(f+g*X(1,k)^2)*X(2,k)+sin(a)*(c5+h*X(1,k)^2)*X(1,k);];
            % fprintf('k=%d p=%.2f w=%.2f w0=%.2f\n',k,X(1,k),X(2,k),w0);
            k = k+1;
            if(k>1000/dt)
                X(2,k) = NaN;
                break;
            end
        end
        R = [R [w0; X(2,k)]];
    end
    figure;
    scatter(R(1,:),R(2,:),'b.');
    hold on;
    plot(R(1,:), R(1,:),'r');
    xlabel('\omega_0');
    ylabel('\omega_{fin}');
    title(['Poincare map for \alpha = ' num2str(a)]);
    hold off;
end

