%clear all 
close all
%clc

barLength = 1.65; % SUPERball length, m
lims = barLength*1;

tspan = 0.05;          % time between plot updates in seconds
delT = 0.001;         % timestep for dynamic sim in seconds
delTUKF  = 0.005;
K = 6e3;              % String stiffness (N/m)
nodalMass = 3*ones(12,1); % Each bar weighs 6kg on SUPERball.
c = 40;             % damping constant, too lazy to figure out units.
F = zeros(12, 3);
stringStiffness = K*ones(24,1);
barStiffness = 10000e3*ones(6,1);
stringDamping = c*ones(24,1);  %string damping vector

addpath('../tensegrityObjects'); % Note forward slash means *nix OS.

% Define SUPERball nodes using Sultan's ordering:
% Initial starting configuration:
l = barLength;
b = l*sqrt(3.0/8.0);
%b = l * 0.9;
%delta = acos(1.0/sqrt(3.0));
delta = 70.0/180.0*pi;
alpha = 60/180*pi;

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

nodes(:,3) = nodes(:,3) - min(nodes(:,3)); % Make minimum node z=0 height.

% Define bars and strings:
bars = [1:2:11; 
        2:2:12];
         %|1 Saddle        6  |7 Vertical       12 |13 Diagonal   18 |19 Boundary     24
strings = [6  9  4  7  2   11  1  5  3  11  7   9   1  3  5  2  4  6  1  5  3  8  12  10;
           9  4  7  2  11  6   6  4  2  8   10  12  11 7  9  8 10 12  5  3  1  12 10  8 ];

% Compute rest lengths for no additional forces:
%stringRestLength = eqManiSVDBRestLengths(alpha, ...
%    delta, barLength, b, K);
stringRestLength=1*ones(24,1); % Placeholder rest lengths.

% Create superball:
barRad = 0.05/2; % Bar diameter for superball is 5cm min.
superBallDynamicsPlot = TensegrityPlot(nodes, strings, bars, barRad,0.005);
superBall = TensegrityStructure(nodes, strings, bars, F, stringStiffness,...
    barStiffness, stringDamping, nodalMass, delT, delTUKF, stringRestLength);


%%%%%% Dynamics Plot %%%%%%%%%%%
% Use a method within TensegrityPlot class to generate a plot of the structure

f = figure('units','normalized', 'outerposition',[0 0 1 1]);
generatePlot(superBallDynamicsPlot, gca);
updatePlot(superBallDynamicsPlot);

%settings to make it pretty
axis equal
view(90, 0); % X-Z view
%view(3)
grid on
light('Position',[0 0 10],'Style','local')
lighting flat
colormap([0.8 0.8 1; 0 1 1])
xlim([-lims lims])
ylim([-lims lims])
zlim(1.6*[-0.01 lims])

drawnow; % Draw and hold initial conditions
%pause(2);

%drawnow; % Hold first frame
%pause(2);

% Create an object calls from the class TensegrityCallbackFunctions which is
% just a function wrapper to hold a vector of persistent values, and some update functions
% It also holds a function for creating sliders
vec = [delta, alpha, b];
calls = TensegrityCallbackFunctions(vec);

%Use the function to create some slider below and attach some callback
%functions to update angle, bed axis, twist, and NR
calls.makeSlider(f,0,0,@(val) updateVal(calls,val,1),[0, pi/2],delta,'delta', 180/pi, '%.0f')
calls.makeSlider(f,475,0,@(val) updateVal(calls,val,2),[pi/6, pi/2],alpha,'alpha', 180/pi, '%.0f')
calls.makeSlider(f,2*475,0,@(val) updateVal(calls,val,3),[0, barLength],b,'baselen/barlen', 1/barLength, '%.3f')

%A custom function needed for each tensegrity structure to update the nodes
%however you see fit, essentially the first time I call this function I set
%some persistent objects/structures (all the items after the ... on the second
%line) this speeds things up since less memory is passed to the function
% each call see the actual function for more details
superBallUpdate(vec, superBall, superBallDynamicsPlot, tspan, barLength, nodalMass(1));

%Create a function handle which only passes the vector of values we will be
%updating
graphicsUpdates = @(vec) superBallUpdate(vec);

t = timer;
t.TimerFcn = @(myTimerObj, thisEvent) timerUpdate(calls, graphicsUpdates);
t.Period = tspan; % Display update time interval
t.ExecutionMode = 'fixedRate';
start(t);
