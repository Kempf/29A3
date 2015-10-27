function [  ] = lcyc( a0, af, da, tol, dt, dw, wtol )
%LCYC Plots limit cycle location for a given range of AoE.
%   a0, af - angle of attack range (0, pi/2)
%   da - angle of attack step (pi/250)
%   tol - omega tolerance (0.01)
%   dt - time step of the solutions (0.01)
%   dw - resolution of the p-map (0.01)
%   wtol - limit cycle omega tolerance (0.001)
    % Initialize stability derivatives 
    tau = 0.5;
    c1 = 0.2;
    c2 = -0.1;
    c3 = 0.1;
    c4  = -1;
    c5 = -0.6;
    c6 = -0.02;
    R = [];
    fprintf('\nSimulating...\t');
    str = '';
    for a = a0:da:af
        % Calculate function values
        f = c1*sin(a)+c2*cos(a);
        g = sin(a)*(c3*(sin(a))^2+(c4+c2)*sin(a)*cos(a)-c1/2*(cos(a))^2);
        h = c6*(sin(a))^2-c5/6*(cos(a))^2;
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
            % Detect limit cycle
            if((abs(w0-X(2,k)) < wtol) && (w0 > tol))
                R = [R [a; X(2,k)]];
            else
                R = [R [a; NaN]];
            end
            % Print progress
            rem = repmat('\b',1,length(str)-1);
            str = [num2str((a-a0)/(af-a0)*100,'%.0f') '%%'];
            fprintf([rem str]);
        end
    end
    fprintf('\n Done!\n');
    close all;
    figure;
    scatter(R(1,:),R(2,:),'b.');
    xlabel('\alpha');
    ylabel('\omega_{fin}');
    title(['Limit cycles for ' num2str(a0) ' \leq \alpha \leq ' num2str(af)]);
end

