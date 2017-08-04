% Plot the projection error of p onto range(A).

npts = 80;

% Connection matrix:
C = [0,0,0,0,0,1,0,0,-1,0,0,0;0,0,0,1,0,0,0,0,-1,0,0,0;0,0,0,1,0,0,-1,0,0,0,0,0;0,1,0,0,0,0,-1,0,0,0,0,0;0,1,0,0,0,0,0,0,0,0,-1,0;0,0,0,0,0,1,0,0,0,0,-1,0;1,0,0,0,0,-1,0,0,0,0,0,0;0,0,0,1,-1,0,0,0,0,0,0,0;0,1,-1,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,1,0,0,-1,0;0,0,0,0,0,0,1,0,0,-1,0,0;0,0,0,0,0,0,0,0,1,0,0,-1;1,0,0,0,0,0,0,0,0,0,-1,0;0,0,1,0,0,0,-1,0,0,0,0,0;0,0,0,0,1,0,0,0,-1,0,0,0;0,1,0,0,0,0,0,-1,0,0,0,0;0,0,0,1,0,0,0,0,0,-1,0,0;0,0,0,0,0,1,0,0,0,0,0,-1;1,0,0,0,-1,0,0,0,0,0,0,0;0,0,1,0,-1,0,0,0,0,0,0,0;1,0,-1,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,1,0,0,0,-1;0,0,0,0,0,0,0,0,0,1,0,-1;0,0,0,0,0,0,0,1,0,-1,0,0;1,-1,0,0,0,0,0,0,0,0,0,0;0,0,1,-1,0,0,0,0,0,0,0,0;0,0,0,0,1,-1,0,0,0,0,0,0;0,0,0,0,0,0,1,-1,0,0,0,0;0,0,0,0,0,0,0,0,1,-1,0,0;0,0,0,0,0,0,0,0,0,0,1,-1];
s = 24; % Number of strings
r = 6;  % Number of rods

% External force vector
p = [zeros(24,1); [3 -1 3 -1 3 -1 -1 -1 -1 -1 -1 -1]'] * 3*9.81; % Each node is 3kg
p = [zeros(24,1); [1 -1 1 -1 1 -1   1 -1  1 -1  1 -1]'] * 3 * 9.81;

% Fixed parameters:
% Define SUPERball nodes using Sultan's ordering:
l = 1;
b = l*sqrt(3.0/8.0);

%for b = linspace(0.1, 1, 3)

alpha_i = linspace(0, 90, npts);
delta_i = linspace(0, 89, npts);
[Alph,Delt] = meshgrid(alpha_i, delta_i);
Z = Alph*0 + NaN;
Z_feas = Alph*0 + NaN;

    for x = 1:npts
        for y = 1:npts
            alpha = alpha_i(x)/180*pi;
            delta = delta_i(y)/180*pi;
            
            if l*sin(delta)*abs(cos(alpha + pi/6)) < b/(2*sqrt(3)) ...
                && sin(alpha + pi/6) < 3*l*sin(delta)/(2*b)
                    ; % We are on eq manifold
            else
                continue; % Point not on eq manifold
            end

            u = sin(delta) * cos(alpha + pi/6.0);
            h = cos(delta)/(2.0*u) * (l*u + sqrt(b*b/3.0 - 3.0*l*l*u*u) - b/sqrt(3.0));
            
            % Handle discontinuity at alpha=pi/3
            if abs(alpha - pi/3.0) < 0.001
                h = l * cos(delta)/2.0;
            end
            
            nodes = [ ...
              % Bar 11, nodes 0, 1
              getNodesForBarWithComAzDec( [l/2.0*sin(delta)*cos(alpha)-b/2.0;
                                l/2.0*sin(delta)*sin(alpha) - b/(2.0*sqrt(3.0));
                                l/2.0*cos(delta)], ...
                                alpha, delta, l);
              % Bar 21, nodes 2, 3
              getNodesForBarWithComAzDec( [l/2.0*sin(delta)*cos(alpha + 4.0*pi/3.0);
                                        b/sqrt(3.0) + l/2.0*sin(delta)*sin(alpha + 4.0*pi/3.0);
                                        l/2.0*cos(delta)], ...
                                        alpha + 4.0*pi/3.0, delta, l);
              % Bar 31, nodes 4, 5
              getNodesForBarWithComAzDec( [b/2.0 + l/2.0*sin(delta)*cos(alpha + 2.0*pi/3.0);
                                        l/2.0*sin(delta)*sin(alpha + 2.0*pi/3.0) - b/(2.0*sqrt(3.0));
                                        l/2.0*cos(delta)], ...
                                        alpha + 2.0*pi/3.0, delta, l);
              % Bar 12, nodes 6, 7
              getNodesForBarWithComAzDec( [l/4*sin(delta)*cos(alpha) + sqrt(3.0)/4.0*l*sin(delta)*sin(alpha) - b/2.0;
                 b/(2.0*sqrt(3.0)) - sqrt(3.0)/4.0*l*sin(delta)*cos(alpha) + l/4.0*sin(delta)*sin(alpha);
                                        3.0/2.0*l*cos(delta) - h], ...
                                        alpha + 2.0*pi/3.0, delta, l);
              % Bar 22, nodes 8, 9
              getNodesForBarWithComAzDec( [b/2.0 - l/2.0*sin(delta)*cos(alpha);
                                        b/(2.0*sqrt(3.0)) - l/2.0*sin(delta)*sin(alpha);
                                        3.0/2.0*l*cos(delta) - h], ...
                                        alpha, delta, l);
              % Bar 32, nodes 10, 11
              getNodesForBarWithComAzDec( [l/4.0*sin(delta)*cos(alpha) - sqrt(3.0)/4.0*l*sin(delta)*sin(alpha);
                       l/4.0*sin(delta)*sin(alpha) + sqrt(3.0)/4.0*l*sin(delta)*cos(alpha) - b/sqrt(3.0);
                                        3.0/2.0*l*cos(delta) - h], ...
                                        alpha + 4.0*pi/3.0, delta, l);            ];
                                    
            % Compute A ,the equilibrium matrix
            A = [transpose(C) * diag(C*nodes(:,1));
                 transpose(C) * diag(C*nodes(:,2));
                 transpose(C) * diag(C*nodes(:,3))];

            %if any(pinv(A)*p < 0)
            %    pinv(A)*p
            %end
            projectionError = transpose(A*pinv(A)*p - p) * (A*pinv(A)*p - p) / norm(p)^2;
            Z(x, y) = sqrt(projectionError);
              
            
%             % Now, solve for force density vector q using quadratic program (QP):
%             Ainv = pinv(A);
%             AinvA = Ainv * A;
% 
%             % Full V.
%             V = (eye(length(AinvA)) - AinvA);
%             
%             % Find an orthonormal basis for the range of V.
%             % NOTE: We don't use orth here, because orth uses precision
%             % eps, which is too small!
%             [Q,R,E] = qr(V);
%             [m , n] = size(R);
%             j=1;
%             i=1;
%             while i<=m
%                 if norm(R(i,:))>10^-12
%                     R_new(j,:)=R(i,:);
%                     j=j+1;
%                 else
%                     i=m;
%                 end
%                 i=i+1;
%             end
%             V=Q(:,1:j-1);
%             
%             c = 0; % Minimum force density for strings.
%             f = 2*transpose(p)*transpose(Ainv(1:s, :))*V(1:s, :);
%             w = quadprog(2*transpose(V(1:s,:))*V(1:s, :), f, -V(1:s, :), zeros(s,1) - c + Ainv(1:s,:)*p);
% 
%             % Now w has been found, we can find q_s:
%             q_s = Ainv(1:s, :) * p + V(1:s, :) * w;
% 
%             % Find q_r, the force densities in the rods, to balance q_s:
%             A_s = A(:, 1:s);
%             A_r = A(:, s+1:end);
%             q_r = pinv(A_r)*(p - A_s*q_s);
% 
%             % We now have q_total
%             q = [q_s; q_r];
%             proj_feas = sqrt(transpose(A*q-p)*(A*q-p) / norm(p)^2);
%             
%             if proj_feas > 2
%                 proj_feas = 2;
%             end
%             Z_feas(x, y) = proj_feas;
        end
    end

surf(Alph, Delt, Z); hold on;
%surf(Alph, Delt, Z_feas);
%contourf(Alph, Delt, Z, 20);
ylabel('alpha');
ylim([30, 90]);
xlabel('delta');
colorbar;

zlabel('error');

pause(1);
drawnow;
%end