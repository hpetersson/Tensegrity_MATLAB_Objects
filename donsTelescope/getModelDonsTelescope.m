function [donsTelescopeModel, nodes, strings, bars] = getModelDonsTelescope()
% This function outputs the model parameters for Don's Telescope
%   Define tensegrity structure

barSpacing = 0.375;
barLength = 1.7;
nodes = [-barSpacing     barLength*0.5  0;
         -barSpacing    -barLength*0.5  0;
          barSpacing     barLength*0.5  0;
          barSpacing    -barLength*0.5  0;
          0             -barSpacing     barLength*0.5;
          0             -barSpacing    -barLength*0.5;
          0              barSpacing     barLength*0.5;
          0              barSpacing    -barLength*0.5;        
          barLength*0.5  0             -barSpacing;
         -barLength*0.5  0             -barSpacing;
          barLength*0.5  0              barSpacing;
         -barLength*0.5  0              barSpacing];
     
HH  = makehgtform('axisrotate',[1 1 0],0.3);
nodes = (HH(1:3,1:3)*nodes')';
nodes(:,3) = nodes(:,3) + 1*barLength;
     
bars = [1:2:11; 
        2:2:12];
strings = [1  1   1  1  2  2  2  2  3  3  3  3  4  4  4  4  5  5  6  6  7  7  8  8;
           7  8  10 12  5  6 10 12  7  8  9 11  5  6  9 11 11 12  9 10 11 12  9 10];

stringRestLength = 0.9*ones(24,1)*norm(nodes(1,:)-nodes(7,:));

K = 1000;
c = 40;
stringStiffness = K*ones(24,1); % String stiffness (N/m)
barStiffness = 100000*ones(6,1); % Bar stiffness (N/m)
stringDamping = c*ones(24,1);  % String damping vector
nodalMass = 1.625*ones(12,1);
delT = 0.001;

donsTelescopeModel = TensegrityStructure(nodes, strings, bars, zeros(12,3), stringStiffness,...
    barStiffness, stringDamping, nodalMass, delT, delT, stringRestLength);
end

