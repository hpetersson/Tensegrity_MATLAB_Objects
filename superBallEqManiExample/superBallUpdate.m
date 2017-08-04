function superBallUpdate(vec, superBall1, superBallDynamicsPlot1,tspan1, ...
    barlength1, nodalmass1)
% Display update function for tensegrity simulation.
global stringTensionDatastore;
global stringTensionCmdDataStore;
global restLenDataStore;
global actualLenCmdLenDataStore;
global HH;

%create some persistent variables for objects and structs
persistent superBall superBallDynamicsPlot tspan barlen baselen alpha_i delta_i nit nodalmass isPacking

if nargin>1
    delta_i = vec(1);
    alpha_i = vec(2);
    baselen = vec(3);
    
    isPacking = 1;
    
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

%%%%%% Update control inputs %%%%%%
delta = delta_i;
alpha = alpha_i;
b = baselen;

% Constant controller
l = barlen;
% Check if point is on eq. manifold:
if l*sin(vec(1))*abs(cos(vec(2) + pi/6)) < vec(3)/(2*sqrt(3)) ...
    && sin(vec(2) + pi/6) < 3*l*sin(vec(1))/(2*vec(3))
    % On manifold!
    delta_i = vec(1);
    alpha_i = vec(2);
    baselen = vec(3);
    delta = vec(1);
    alpha = vec(2);
    b = vec(3);
else
    % Not on manifold!
    fprintf('Not on mani!\n');
end

% Control alpha(t), delta(t):
nit = nit + 1;
if isPacking
    time = nit * tspan - 1;
else
    time = 15-nit*tspan;
end
if time < 0
   time = 0;
end

if time > 9
    time=9;
end

% if isPacking && time > 15
%     nit = 0;
%     isPacking = 0;
% end

% Commanded alpha(t), delta(t):
Td=10;
%Cylindrical packing:
delta = max(9/180*pi, delta_i + time/Td*(9/180*pi-delta_i));
alpha = min(70/180*pi, alpha_i + time/Td*(80/180*pi-alpha_i)); % should be 70
b = max(0.1*barlen, baselen + (time)/Td*(0.1*barlen - baselen)); % should be 0.1
if b > baselen
    b = baselen;
end



% Deployment from Star Packing:
% delta = max(55/180*pi, delta_i + time/Td*(55/180*pi-delta_i))
% alpha = 60/180*pi;
% b = baselen;

% Star Packing:
% delta = min(88/180*pi, delta_i + time/Td*(88/180*pi-delta_i));
% alpha = 60/180*pi;
% b = baselen;

% delta = 55/180*pi;
% alpha = 60/180*pi;
% b = baselen;

%alpha/pi*180
%delta/pi*180
%b
%c1 = b./(2*sqrt(3)*l*sin(delta));
%acosc = acos(c1);
% %acosc(imag(acosc)~=0) = NaN;
% %plot(delta_i, (-acosc-pi/6 + pi)/pi*180); hold on;
%alpha = -acosc-pi/6 + pi - 0.1; % Alpha on edge of eq manifold


%%%%% Compute commanded rest lengths: %%%%%

% Define correction forces:
C = superBall.C;
s = 24; % Number of strings
r = 6;  % Number of rods

nodes = SVDB2nodes(alpha, delta, b, l);
% Rotate so desired face is on ground:
nodes = (HH(1:3,1:3)*nodes')';

% Compute A ,the equilibrium matrix
A = [transpose(C) * diag(C*nodes(:,1));
     transpose(C) * diag(C*nodes(:,2));
     transpose(C) * diag(C*nodes(:,3))];
 
% Now, solve for force density vector q using quadratic program (QP):
Ainv = pinv(A, 1e-3);
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

%c = 40; % Minimum force density for strings.
c = 10;
%c = 0;
% Desired external forces
%c = 1;

% p for "rolling motion":
% nodest = nodes;
% nodest(:,3) = nodest(:,3) - mean(nodest(:,3));
% HH  = makehgtform('axisrotate',[1 0 0], max(0, time/2 - 2));
% nodest = (HH(1:3,1:3)*nodest')';
% z = nodest(:,3);
% [~,min1] = min(z); z(min1) = 1e9;
% [~,min2] = min(z); z(min2) = 1e9;
% [~,min3] = min(z); z(min3) = 1e9;
% touchingGround = z>1e8;
% pz = -1*ones(12,1);
% pz(touchingGround) = pz(touchingGround) + 12/sum(double(touchingGround));
% p = [zeros(24,1); pz] * nodalmass * 9.81;
% HH2  = makehgtform('axisrotate',[1 0 0], 1);
% pt = (HH2(1:3,1:3)*[p(1:12), p(13:24), p(25:end)]')';
% p = [pt(:,1); pt(:,2); pt(:,3)]*2;

% For (0,4,2)+1 face on ground:
p = [zeros(24,1); [3 -1 3 -1 3 -1  -1 -1 -1 -1 -1 -1]'] * nodalmass * 9.81;
% For (10,5,0)+1 face on ground:
%p = [zeros(24,1); [3 -1 -1 -1 -1 3  -1 -1 -1 -1 3 -1]'] * nodalmass * 9.81;

%p = zeros(12*3,1);
f = 2*transpose(p)*transpose(Ainv(1:s, :))*V(1:s, :);
options = optimoptions('quadprog','Algorithm',  'interior-point-convex','Display','off');
w = quadprog(2*transpose(V(1:s,:))*V(1:s, :), f, -V(1:s, :), zeros(s,1) - c + Ainv(1:s,:)*p, ...
    [],[], [],[], [], options);

% Now w has been found, we can find q_s:
q_s = Ainv(1:s, :) * p + V(1:s, :) * w;
%q_s_particular = Ainv(1:s, :) * p;
%q_s_null = V(1:s, :) * w;
% Limit particular solution if it's too big!
% if norm(q_s_particular)/norm(q_s_null) > 0.3
%    q_s_particular = q_s_particular / norm(q_s_particular) * norm(q_s_null)*0.3;
% end
% if max(q_s_particular) > 170
%     q_s_particular = q_s_particular/max(q_s_particular)*170;
% end
% q_s = q_s_particular + q_s_null/max(q_s_null)*(200 - max(q_s_particular));
%if max(q_s_particular) > 200
%    q_s_particular = q_s_particular/max(q_s_particular)

% Find q_r, the force densities in the rods, to balance q_s:
A_s = A(:, 1:s);
A_r = A(:, s+1:end);
q_r = pinv(A_r)*(p - A_s*q_s);

% We now have q_total
q = [q_s; q_r];

% Check we have indeed found an equilibrium solution:
error = norm(p - A*q)/norm(p);

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
           ones(6, 1) * B;
           ones(6, 1) * barlen];
% xlens = C*nodes(:,1);
% ylens = C*nodes(:,2);
% zlens = C*nodes(:,3);
% lengths = sum([xlens ylens zlens].^2, 2).^0.5;
% qmin = ones(length(q),1)*10 ./ lengths; % min/max q limits for strings
% qmax = ones(length(q),1)*200 ./ lengths;
% b = 6;
% qmin(end-b+1:end) = -1e6; % No limits on q in bars
% qmax(end-b+1:end) = 1e6;
% q = quadprog(2*(A')*A + 0*eye(length(q))*norm(A)/1000,... % H
%          -2*p'*A,... % f
%          [],... % Inequality constraint matrix
%          [],... % Inequality constraint vector
%          [],... % Equality constraint matrix
%          [],... % Equality constraint vector
%          qmin,... % q lower bound
%          qmax,... % q upper bound
%          [],... % x0
%          options);
% q_s = q(1:s);

T_cmd = q_s .* lengths(1:s); % Only take string force densities
T_cmd(T_cmd>200)=200; % Clip tensions
%if max(T_cmd) > 200
%    T_cmd = T_cmd/max(T_cmd)*200; % Limit tensions to 200N
%end
restlens = lengths(1:s) - T_cmd ./ superBall.simStruct.stringStiffness; % New rest lengths
%restlens = restlens + 0.1/(time+0.1); % "Ease-in"
%restlens=[0.548476888595409;1.70512319073480;0.577041899562893;1.71007669733963;0.519994217307759;1.66317996617458;1.27361609504908;1.29788434153900;1.29187713658366;0.347298492599110;0.344538497008120;0.347637055448560;1.12175240752246;1.12760356730137;1.08598285145707;0.430445297174641;0.342371995645333;0.456788473629911;0.275963850329595;0.252318995886744;0.234837866982215;1.27598866123661;1.28924970275127;1.38880430178505];
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

