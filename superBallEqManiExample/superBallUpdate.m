function superBallUpdate(superBall1,superBallDynamicsPlot1,tspan1, ...
    barlength, baselength, alphai, deltai)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%create some persistent variables for objects and structs
persistent superBall superBallDynamicsPlot tspan barlen baselen alpha_i delta_i nit

if nargin>1
    superBall = superBall1;
    superBallDynamicsPlot = superBallDynamicsPlot1;
    tspan = tspan1;
    barlen = barlength;
    baselen = baselength;
    alpha_i = alphai;
    delta_i = deltai;
    nit = 0;
end

%%%%%%%%%%%% Update dynamics $%%%%%%%%%%%%%%%%%

dynamicsUpdate(superBall, tspan);
actualNodes = superBall.ySim(1:end/2,:);
%barVec = actualNodes(bars(1,:),:) - actualNodes(bars(2,:),:);
%barNorm = sqrt(barVec(:,1).^2 + barVec(:,2).^2 + barVec(:,3).^2);
%barAngleFromVert = acos(barVec(:,3:3:end)./barNorm);
superBallDynamicsPlot.nodePoints = actualNodes;
updatePlot(superBallDynamicsPlot);
drawnow  %plot it up

nit = nit + 1;
time = nit * tspan;

% Commanded alpha(t), delta(t):
delta = (55.0-delta_i/pi*180)*(time-2.0)/5.0 + delta_i/pi*180;
if delta < 55.0
    delta = 55.0;
end
if delta > 70.0 
    delta = 70.0;
end

delta = delta/180.0 * pi; % rad

%delta = delta_i; % Constant controller

% Update rest length commands:
%superBall.simStruct.stringRestLengths = eqManiSVDBRestLengths(alpha_i, ...
%    delta, barlen, baselen, superBall.simStruct.stringStiffness(1));
end

