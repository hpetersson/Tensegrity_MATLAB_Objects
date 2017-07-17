%clear all 
close all
%clc

barLength = 1.65; % SUPERball length
lims = barLength*1;

tspan = 0.05;          % time between plot updates in seconds
delT = 0.001;         % timestep for dynamic sim in seconds
delTUKF  = 0.005;
K = 1e3;              % String stiffness (N/m)
nodalMass = 3*ones(12,1); % Each bar weighs 6kg on SUPERball.
c = 40;             % damping constant, too lazy to figure out units.
F = zeros(12, 3);
stringStiffness = K*ones(24,1);
barStiffness = 100e3*ones(6,1);
stringDamping = c*ones(24,1);  %string damping vector

addpath('../tensegrityObjects'); % Note forward slash means *nix OS.

% Define SUPERball nodes using Sultan's ordering:
% Initial starting configuration:
l = barLength;
b = l*sqrt(3.0/8.0);
%b = l * 0.9;
%delta = acos(1.0/sqrt(3.0));
delta = 85.0/180.0*pi;
alpha = 68/180*pi;

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
  getNodesForBarWithComAzDec( [l/4.0*sin(delta)*cos(alpha) - sqrt(3.0)/4.0*l*sin(delta)*sin(alpha),
           l/4.0*sin(delta)*sin(alpha) + sqrt(3.0)/4.0*l*sin(delta)*cos(alpha) - b/sqrt(3.0),
                            3.0/2.0*l*cos(delta) - h], ...
                            alpha + 4.0*pi/3.0, delta, l);            ];
                            
%HH  = makehgtform('axisrotate',[1 1 0],0.3);
%     nodes = (HH(1:3,1:3)*nodes')';

nodes(:,3) = nodes(:,3) - min(nodes(:,3)) ; % Make minimum node z=0 height.

% Define bars and strings:
bars = [1:2:11; 
        2:2:12];
         %|1 Saddle        6  |7 Vertical       12 |13 Diagonal   18 |19 Boundary     24
strings = [6  9  4  7  2   11  1  5  3  11  7   9   1  3  5  2  4  6  1  5  3  8  12  10;
           9  4  7  2  11  6   6  4  2  8   10  12  11 7  9  8 10 12  5  3  1  12 10  8 ];

% Compute rest lengths for no additional forces:
%stringRestLength = eqManiSVDBRestLengths(alpha, ...
%    delta, barLength, b, K);
stringRestLength=1*ones(24,1);

barRad = 0.05/2; % Bar diameter for superball is 5cm min.
superBallDynamicsPlot = TensegrityPlot(nodes, strings, bars, barRad,0.005);
superBall = TensegrityStructure(nodes, strings, bars, F, stringStiffness,...
    barStiffness, stringDamping, nodalMass, delT, delTUKF, stringRestLength);


% Define correction forces:
C = superBall.C;
s = 24; % Number of strings
r = 6;  % Number of rods

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

c = 200; % Minimum force density for strings.
% Desired external forces
p = [zeros(24,1); [3 -1 3 -1 3 -1 -1 -1 -1 -1 -1 -1]'] * nodalMass(1)*9.81;
%p = zeros(12*3,1);
f = 2*transpose(p)*transpose(Ainv(1:s, :))*V(1:s, :);
w = quadprog(2*transpose(V(1:s,:))*V(1:s, :), f, -V(1:s, :), zeros(s,1) - c + Ainv(1:s,:)*p);
w = 100;
% Now w has been found, we can find q_s:
q_s = Ainv(1:s, :) * p + V(1:s, :) * w

% Find q_r, the force densities in the rods, to balance q_s:
A_s = A(:, 1:s);
A_r = A(:, s+1:end);
q_r = pinv(A_r)*(p - A_s*q_s);

% We now have q_total
q = [q_s; q_r];

% Check we have indeed found an equilibrium solution:
error = norm(p - A*q)/norm(p)

S = sqrt(h*h + b*b/3.0 + l*l*sin(delta)*sin(delta) ...
           - 2.0/sqrt(3.0)*l*b*sin(delta)*cos(alpha - pi/6.0) );
V = sqrt(b*b + l*l - 2.0*l*b*sin(delta)*sin(alpha + pi/6.0));
D = sqrt(h*h + b*b/3.0 + l*l - 2.0/sqrt(3.0)*l*b*sin(delta)*sin(alpha) - 2.0*l*h*cos(delta));
B = b;

% Handle discontinuity at alpha=pi/3
if abs(alpha - pi/3.0) < 0.001
    h = l * cos(delta)/2.0;
    S = sqrt(h*h + b*b/3.0 + l*l*sin(delta)*sin(delta) ...
           - 2.0/sqrt(3.0)*l*b*sin(delta)*cos(alpha - pi/6.0) );
    V = sqrt(b*b + l*l - 2.0*l*b*sin(delta)*sin(alpha + pi/6.0));
    D = sqrt(h*h + b*b/3.0 + l*l - 2.0/sqrt(3.0)*l*b*sin(delta)*sin(alpha) - 2.0*l*h*cos(delta));
end

lengths = [ones(6, 1) * S;
           ones(6, 1) * V;
           ones(6, 1) * D;
           ones(6, 1) * B ];
T_str = q(1:24) .* lengths % Only take string force densities
%T_str = T_str/norm(T_str)*200
restlens = lengths - T_str / K;

% Update rest lengths
%restlens = (restlens-lengths)*0.9 + eqManiSVDBRestLengths(alpha, delta, barLength, b, K)
superBall.simStruct.stringRestLengths = restlens;

%%%%%% Dynamics Plot %%%%%%%%%%%
% Use a method within TensegrityPlot class to generate a plot of the structure

f = figure('units','normalized', 'outerposition',[0 0 1 1]);
generatePlot(superBallDynamicsPlot, gca);
updatePlot(superBallDynamicsPlot);

%settings to make it pretty
axis equal
view(90, 0)
%view(3)
grid on
light('Position',[0 0 10],'Style','local')
lighting flat
colormap([0.8 0.8 1; 0 1 1])
xlim([-lims lims])
ylim([-lims lims])
zlim(1.6*[-0.01 lims])

drawnow; % Draw and hold initial conditions
pause(2);
superBallUpdate(superBall, superBallDynamicsPlot, tspan, barLength, b, ...
    alpha, delta)

drawnow; % Hold first frame
%pause(2);


for i = 1:1000e3
    superBallUpdate
    %pause(0.1);
  %  MM(i) = getframe(f);
end
%filename = 'quickAnimation.avi';
%writerObj = VideoWriter(filename);
%writerObj.FrameRate = 20;
%open(writerObj);
%writeVideo(writerObj,MM);
%close(writerObj);
% 
% t = timer;
% t.TimerFcn = @(myTimerObj, thisEvent) superBallUpdate;
% t.Period = tspan;
% t.ExecutionMode = 'fixedRate';
% start(t);

% % 

