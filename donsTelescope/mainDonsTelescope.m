% Test to set up Don Peckhams tensegrity telescope

clear all;
close all;
addpath('../tensegrityObjects')

%% Set up model 
% Define tensegrity structure

% Three types of bars: horizontal (usually in a square configuration);
% vertical (holding up the secondary mirror); and circular (at bottom and
% top)

barSpacing = 0.375;
barLengthVertical = 1.5;
barDistVert = 0.3;
barLengthHorizontal = 0.4;

nodes = [ barDistVert   0               0;                  
          barDistVert   barDistVert     0;
          0             barDistVert     0;
          0             0               0;          
          barDistVert   0               barLengthVertical;
          barDistVert   barDistVert     barLengthVertical;
          0             barDistVert     barLengthVertical;
          0             0               barLengthVertical;
          0.35          -0.05           barLengthVertical/2;
          0.35          0.35            barLengthVertical/2;
          -0.05         0.35            barLengthVertical/2;
          -0.05         -0.05           barLengthVertical/2
                                        ]
     
HH  = makehgtform('axisrotate',[1 1 0],0.0);
nodes = (HH(1:3,1:3)*nodes')';

     
bars = [1 2 3 4 9 10 11 12 1 2 3 4 5 6 7 8; 
        5 6 7 8 10 11 12 9 2 3 4 1 6 7 8 5];
strings = [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8;
           10 12 9 11 10 12 9 11 10 12 9 11 10 12 9 11];

stringRestLength = 0.9*ones(24,1)*norm(nodes(1,:)-nodes(7,:));

K = 1000;
c = 40;
stringStiffness = K*ones(24,1); % String stiffness (N/m)
barStiffness = 100000*ones(6,1); % Bar stiffness (N/m)
stringDamping = c*ones(24,1);  % String damping vector
nodalMass = 1.625*ones(12,1);
delT = 0.001;

donsTelescopeModel = TensegrityStructure(nodes, strings, bars, zeros(numel(nodes(:,1)),3), stringStiffness,...
    barStiffness, stringDamping, nodalMass, delT, delT, stringRestLength);

%% Create dynamics display
bar_radius = 0.025; % meters
string_radius = 0.005;
superBallDynamicsPlot = TensegrityPlot(nodes, strings, bars, bar_radius, string_radius);
f = figure('units','normalized','outerposition',[0 0 1 1]);
% use a method within TensegrityPlot class to generate a plot of the
% structure
generatePlot(superBallDynamicsPlot,gca);
updatePlot(superBallDynamicsPlot);

%settings to make it pretty
axis equal
view(3)
grid on
light('Position',[0 0 10],'Style','local')
lighting flat
colormap([0.8 0.8 1; 0 1 1]);
lims = 1.2*barLengthVertical;
xlim([-lims lims])
ylim([-lims lims])
zlim(1.6*[-0.01 lims])

%% Run dynamics
displayTimespan = 0.05; % 20fps. Increase display time interval if system can't keep up.
myDynamicsUpdate(donsTelescopeModel, superBallDynamicsPlot, displayTimespan);

for i = 1:200
    myDynamicsUpdate();
end