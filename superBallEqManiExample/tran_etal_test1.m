% Tran et al. formfinding method
clear all
clc

% Connectivity matrix:
C = [0,0,0,0,0,1,0,0,-1,0,0,0;0,0,0,1,0,0,0,0,-1,0,0,0;0,0,0,1,0,0,-1,0,0,0,0,0;0,1,0,0,0,0,-1,0,0,0,0,0;0,1,0,0,0,0,0,0,0,0,-1,0;0,0,0,0,0,1,0,0,0,0,-1,0;1,0,0,0,0,-1,0,0,0,0,0,0;0,0,0,1,-1,0,0,0,0,0,0,0;0,1,-1,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,1,0,0,-1,0;0,0,0,0,0,0,1,0,0,-1,0,0;0,0,0,0,0,0,0,0,1,0,0,-1;1,0,0,0,0,0,0,0,0,0,-1,0;0,0,1,0,0,0,-1,0,0,0,0,0;0,0,0,0,1,0,0,0,-1,0,0,0;0,1,0,0,0,0,0,-1,0,0,0,0;0,0,0,1,0,0,0,0,0,-1,0,0;0,0,0,0,0,1,0,0,0,0,0,-1;1,0,0,0,-1,0,0,0,0,0,0,0;0,0,1,0,-1,0,0,0,0,0,0,0;1,0,-1,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,1,0,0,0,-1;0,0,0,0,0,0,0,0,0,1,0,-1;0,0,0,0,0,0,0,1,0,-1,0,0;1,-1,0,0,0,0,0,0,0,0,0,0;0,0,1,-1,0,0,0,0,0,0,0,0;0,0,0,0,1,-1,0,0,0,0,0,0;0,0,0,0,0,0,1,-1,0,0,0,0;0,0,0,0,0,0,0,0,1,-1,0,0;0,0,0,0,0,0,0,0,0,0,1,-1];
s = 24;
b = 6;
nn = 12; % Number of nodes

q = [ones(s,1); -ones(b,1)];

for niter = 1:4
    D = C'*diag(q)*C;
    [phi, lambda] = eig(D);

    basis = [phi(:,1), phi(:,2), phi(:,3), phi(:,4)];
    Ld = C*basis;
    Ldn = sum(Ld.^2,1);
    nbasis = basis(:, Ldn > 1e-3);
    %if size(nbasis, 2) ~= 3
   % 	nbasis = basis(:, 1:3);
    %end

    n = nbasis(:);
    A = [transpose(C) * diag(C*n(1:nn));
         transpose(C) * diag(C*n(nn+1:2*nn));
         transpose(C) * diag(C*n(2*nn+1:end))];

    [U, V, W] = svd(A);
    q = W(:, end)

    error = norm(A*q)
end

xlens = C*n(1:nn);
ylens = C*n(nn+1:2*nn);
zlens = C*n(2*nn+1:end);
lens = sum([xlens ylens zlens].^2, 2).^0.5

% Display tensegrity:
nodes = [n(1:nn) n(nn+1:2*nn) n(2*nn+1:end)];
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
stringRestLength=2*ones(24,1); % Placeholder rest lengths.

% Create superball:
barRad = 0.02/2; % Bar diameter for superball is 5cm min.
stringRad = 0.002/2;
superBallplot = TensegrityPlot(nodes, strings, bars, barRad, stringRad);
f = figure('units','normalized', 'outerposition',[0 0 1 1]);
generatePlot(superBallplot, gca);
updatePlot(superBallplot);
%settings to make it pretty
axis equal
view(90, 0); % X-Z view
%view(3)
grid on
light('Position',[0 0 10],'Style','local')
lighting flat
colormap([0.8 0.8 1; 0 1 1])
lims = 1;
xlim([-lims lims])
ylim([-lims lims])
zlim(1.6*[-0.01 lims])