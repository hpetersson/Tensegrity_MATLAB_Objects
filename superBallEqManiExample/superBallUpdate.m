function superBallUpdate(vec, superBall1, superBallDynamicsPlot1,tspan1, ...
    barlength1, nodalmass1)
% Display update function for tensegrity simulation.

%create some persistent variables for objects and structs
persistent superBall superBallDynamicsPlot tspan barlen baselen alpha_i delta_i nit nodalmass

if nargin>1
    delta_i = vec(1);
    alpha_i = vec(2);
    baselen = vec(3);
    
    superBall = superBall1;
    superBallDynamicsPlot = superBallDynamicsPlot1;
    tspan = tspan1;
    barlen = barlength1;
    nodalmass = nodalmass1;
    
    nit = 0;
end

%%%%%% Update control inputs %%%%%%
% Control alpha(t), delta(t):
% nit = nit + 1;
% time = nit * tspan;
% 
% % Commanded alpha(t), delta(t):
% delta = (55.0-delta_i/pi*180)*(time-2.0)/5.0 + delta_i/pi*180;
% if delta < 55.0
%     delta = 55.0;
% end
% if delta > 70.0 
%     delta = 70.0;
% end

delta = delta_i;
alpha = alpha_i;
b = baselen;

% Constant controller
l = barlen;
% Check if point is on eq. manifold:
%if l*sin(vec(1))*abs(cos(vec(2) + pi/6)) < vec(3)/(2*sqrt(3)) ...
%    && sin(vec(2) + pi/6) < 3*l*sin(vec(1))/(2*vec(3))
    % On manifold!
    delta = vec(1);
    alpha = vec(2);
    b = vec(3)
%else
    % Not on manifold!
%end

%%%%% Compute commanded rest lengths: %%%%%

% Define correction forces:
C = superBall.C;
s = 24; % Number of strings
r = 6;  % Number of rods

nodes = superBall.nodePoints;

% Compute A ,the equilibrium matrix
A = [transpose(C) * diag(C*nodes(:,1));
     transpose(C) * diag(C*nodes(:,2));
     transpose(C) * diag(C*nodes(:,3))];
 
% Now, solve for force density vector q using quadratic program (QP):
Ainv = pinv(A);
AinvA = Ainv * A;

% Full V.
V = (eye(length(AinvA)) - AinvA);
% Find an orthonormal basis for the range of V.
% NOTE: We don't use orth here, because orth uses precision
% eps, which is too small!
[Q,R,E] = qr(V);
[m , n] = size(R);
j=1;
i=1;
while i<=m
    if norm(R(i,:))>10^-12
        R_new(j,:)=R(i,:);
        j=j+1;
    else
        i=m;
    end
    i=i+1;
end
V=Q(:,1:j-1);

% Use V as the first s rows, for only string force densities.

%c = 200*200/315.9879; % Minimum force density for strings.
c = 1000;
% Desired external forces
%p = [zeros(24,1); [3 -1 3 -1 3 -1 -1 -1 -1 -1 -1 -1]'] * nodalmass * 9.81;
p = zeros(12*3,1);
f = 2*transpose(p)*transpose(Ainv(1:s, :))*V(1:s, :);
options = optimoptions('quadprog','Algorithm',  'interior-point-convex','Display','off');
w = quadprog(2*transpose(V(1:s,:))*V(1:s, :), f, -V(1:s, :), zeros(s,1) - c + Ainv(1:s,:)*p);

% Now w has been found, we can find q_s:
q_s = Ainv(1:s, :) * p + V(1:s, :) * w;

% Find q_r, the force densities in the rods, to balance q_s:
% A_s = A(:, 1:s);
% A_r = A(:, s+1:end);
% q_r = pinv(A_r)*(p - A_s*q_s);
% 
% % We now have q_total
% q = [q_s; q_r];

% Check we have indeed found an equilibrium solution:
%error = norm(p - A*q)/norm(p)

% Compute actual lengths of strings:
u = sin(delta) * cos(alpha + pi/6.0);
h = cos(delta)/(2.0*u) * (l*u + sqrt(b*b/3.0 - 3.0*l*l*u*u) - b/sqrt(3.0));

% Handle discontinuity at alpha=pi/3
if abs(alpha - pi/3.0) < 0.001
    h = l * cos(delta)/2.0;
end

S = sqrt(h*h + b*b/3.0 + l*l*sin(delta)*sin(delta) ...
           - 2.0/sqrt(3.0)*l*b*sin(delta)*cos(alpha - pi/6.0) );
VV = sqrt(b*b + l*l - 2.0*l*b*sin(delta)*sin(alpha + pi/6.0));
D = sqrt(h*h + b*b/3.0 + l*l - 2.0/sqrt(3.0)*l*b*sin(delta)*sin(alpha) - 2.0*l*h*cos(delta));
B = b;

% Handle discontinuity at alpha=pi/3
if abs(alpha - pi/3.0) < 0.001
    S = sqrt(h*h + b*b/3.0 + l*l*sin(delta)*sin(delta) ...
           - 2.0/sqrt(3.0)*l*b*sin(delta)*cos(alpha - pi/6.0) );
    VV = sqrt(b*b + l*l - 2.0*l*b*sin(delta)*sin(alpha + pi/6.0));
    D = sqrt(h*h + b*b/3.0 + l*l - 2.0/sqrt(3.0)*l*b*sin(delta)*sin(alpha) - 2.0*l*h*cos(delta));
end

lengths = [ones(6, 1) * S;
           ones(6, 1) * VV;
           ones(6, 1) * D;
           ones(6, 1) * B ];
T_str = q_s .* lengths; % Only take string force densities
%T_str(T_str>300)=300; % Clip
T_str
%T_str = T_str./max(T_str)*200 % Limit tensions to 200N
restlens = lengths - T_str ./ superBall.simStruct.stringStiffness % New rest lengths
superBall.simStruct.stringRestLengths = restlens;


%%%%%%%%%%%% Update dynamics $%%%%%%%%%%%%%%%%%
dynamicsUpdate(superBall, tspan);
%%% Plot %%%
actualNodes = superBall.ySim(1:end/2,:);
superBallDynamicsPlot.nodePoints = actualNodes;
updatePlot(superBallDynamicsPlot);
drawnow  %plot it up
end

