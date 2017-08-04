function superBallUpdate(vec, superBall1, superBallDynamicsPlot1,tspan1, ...
    barlength1, nodalmass1)
% Display update function for tensegrity simulation.
global stringTensionDatastore;
global stringTensionCmdDataStore;
global restLenDataStore;
global actualLenCmdLenDataStore;

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
    
    % Clear datastores
    stringTensionDatastore = [];
    stringTensionCmdDataStore = [];
    restLenDataStore = [];
    actualLenCmdLenDataStore = [];
    
    nit = 0;
end

%%%%% Compute commanded rest lengths: %%%%%

% Define correction forces:
C = superBall.C;
s = 24; % Number of strings
r = 6;  % Number of rods

alpha = 60/180*pi;
delta = 55/180*pi;
b = sqrt(3/8)*barlen;
nodes = SVDB2nodes(alpha, delta, b, barlen);

% Compute A, the equilibrium matrix
A = [transpose(C) * diag(C*nodes(:,1));
     transpose(C) * diag(C*nodes(:,2));
     transpose(C) * diag(C*nodes(:,3))];
 
nodesOnGroundL = logical([1 0 1 0 1 0 0 0 0 0 0 0]');
idxx = 1:12;
nodesOnGround = idxx(nodesOnGroundL);
nodesOffGround = idxx(~nodesOnGroundL);
% Make new A for nodes on, off ground: (6 x 30).
reduA = [sum(A(nodesOffGround,:), 1); % x, y, z, force for nodes off ground
        sum(A(nodesOffGround+12,:), 1);
        sum(A(nodesOffGround+24,:), 1);
        sum(A(nodesOnGround,:), 1); % x, y, z, force for nodes on ground
        sum(A(nodesOnGround+12,:), 1);
        sum(A(nodesOnGround+24,:), 1) ];

% Now, solve for q that gives desired F_ext_(on, off) ground:
F_ext_des = [1 0 0, -1 0 0]';
 
% Now, solve for force density vector q using quadratic program (QP):
rAinv = pinv(reduA, 1e-3);
Ainv = pinv(A, 1e-3);
rAinvA = rAinv * reduA;

% Full V.
V = (eye(length(rAinvA)) - rAinvA);
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

% Minimum force density for strings.
c = 10;

% For (0,4,2)+1 face on ground:
%p = [zeros(24,1); [3 -1 3 -1 3 -1  -1 -1 -1 -1 -1 -1]'] * nodalmass * 9.81;
% For (10,5,0)+1 face on ground:
%p = [zeros(24,1); [3 -1 -1 -1 -1 3  -1 -1 -1 -1 3 -1]'] * nodalmass * 9.81;
p = F_ext_des;

%p = zeros(12*3,1);
f = 2*transpose(p)*transpose(rAinv(1:s, :))*V(1:s, :);
options = optimoptions('quadprog','Algorithm',  'interior-point-convex','Display','off');
w = quadprog(2*transpose(V(1:s,:))*V(1:s, :), f, -V(1:s, :), zeros(s,1) - c + rAinv(1:s,:)*p, ...
    [],[], [],[], [], options);

cost = w'*2*transpose(V(1:s,:))*V(1:s, :)*w + f*w

% Now w has been found, we can find q_s:
%q_s = Ainv(1:s, :) * p + V(1:s, :) * w;

q = rAinv * p + V * w;
q_s = q(1:s);

% Find q_r, the force densities in the rods, to balance q_s:
A_s = reduA(:, 1:s);
A_r = reduA(:, s+1:end);
q_r = pinv(A_r)*(p - A_s*q_s);

% We now have q_total
q = [q_s; q_r];

% Check we have indeed found an equilibrium solution:
%error = norm(p - A*q)/norm(p);

xlens = C*nodes(:,1);
ylens = C*nodes(:,2);
zlens = C*nodes(:,3);
lengths = sum([xlens ylens zlens].^2, 2).^0.5;

T_cmd = q_s .* lengths(1:s); % Only take string force densities
%T_cmd(T_cmd>200)=200; % Clip tensions
%if max(T_cmd) > 200
%    T_cmd = T_cmd/max(T_cmd)*200; % Limit tensions to 200N
%end
restlens = lengths(1:s) - T_cmd ./ superBall.simStruct.stringStiffness; % New rest lengths
%restlens = restlens + 0.1/(time+0.1); % "Ease-in"
superBall.simStruct.stringRestLengths = restlens;


%%%%%%%%%%%% Update dynamics $%%%%%%%%%%%%%%%%%
dynamicsUpdate(superBall, tspan);
%%% Plot %%%
actualNodes = superBall.ySim(1:end/2,:);
superBallDynamicsPlot.nodePoints = actualNodes;
T_actual = superBall.stringTensions;
updatePlot(superBallDynamicsPlot);
drawnow  %plot it up

% Save Data:
stringTensionDatastore = [stringTensionDatastore; T_actual'];
stringTensionCmdDataStore = [stringTensionCmdDataStore; T_cmd'];
restLenDataStore = [restLenDataStore; restlens'];
actualLenCmdLenDataStore = [actualLenCmdLenDataStore; lengths(1:s)'];
end

