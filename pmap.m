function [  ] = pmap( a )
%SOLN Solves the wing rock ODE and plots a phase map.
%   a - angle of attack (rad)
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
    R = [];
    for w0 = 0:dt:1
        k = 1;
        while ((X(1,k) < pi) && (X(1,k) > -pi))
            X(:,k+1) = X(:,k)+dt*[X(2,k);tau*(f+g*X(1,k)^2)*X(2,k)+sin(a)*(c5+h*X(1,k)^2)*X(1,k);];      
            k = k+1;
        end
        R = [R [w0; mu*X(2,k)]];
    end
    close all;
    figure;
    axis equal manual;
    scatter(R(1,:),R(2,:),'b.');
    hold on;
    plot(R(1,:), R(1,:),'r');
    xlabel('\omega_0');
    ylabel('\omega_{fin}');
    title(['Poincaré map for \alpha = ' num2str(a)]);
    hold off;
end

