% Stiffness and Stability of SUPERball
clear all;
clc

l = 1.65;
b = sqrt(3/8)*l;
alpha = 60/180*pi;
delta = 70/180*pi;
n = SVDB2nodes(alpha, delta, b, l);
n = n(:); % Make column vector
n = n + (rand(length(n), 1)-0.5)*l/5;
n0 = n;

% Connectivity matrix:
C = [0,0,0,0,0,1,0,0,-1,0,0,0;0,0,0,1,0,0,0,0,-1,0,0,0;0,0,0,1,0,0,-1,0,0,0,0,0;0,1,0,0,0,0,-1,0,0,0,0,0;0,1,0,0,0,0,0,0,0,0,-1,0;0,0,0,0,0,1,0,0,0,0,-1,0;1,0,0,0,0,-1,0,0,0,0,0,0;0,0,0,1,-1,0,0,0,0,0,0,0;0,1,-1,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,1,0,0,-1,0;0,0,0,0,0,0,1,0,0,-1,0,0;0,0,0,0,0,0,0,0,1,0,0,-1;1,0,0,0,0,0,0,0,0,0,-1,0;0,0,1,0,0,0,-1,0,0,0,0,0;0,0,0,0,1,0,0,0,-1,0,0,0;0,1,0,0,0,0,0,-1,0,0,0,0;0,0,0,1,0,0,0,0,0,-1,0,0;0,0,0,0,0,1,0,0,0,0,0,-1;1,0,0,0,-1,0,0,0,0,0,0,0;0,0,1,0,-1,0,0,0,0,0,0,0;1,0,-1,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,1,0,0,0,-1;0,0,0,0,0,0,0,0,0,1,0,-1;0,0,0,0,0,0,0,1,0,-1,0,0;1,-1,0,0,0,0,0,0,0,0,0,0;0,0,1,-1,0,0,0,0,0,0,0,0;0,0,0,0,1,-1,0,0,0,0,0,0;0,0,0,0,0,0,1,-1,0,0,0,0;0,0,0,0,0,0,0,0,1,-1,0,0;0,0,0,0,0,0,0,0,0,0,1,-1];
s = 24;
b = 6;
nn = 12; % Number of nodes
q0 = [ones(s,1); -ones(b,1)];

singvals = [];

for niter = 1:10
    % Compute A ,the equilibrium matrix
    A = [transpose(C) * diag(C*n(1:nn));
         transpose(C) * diag(C*n(nn+1:2*nn));
         transpose(C) * diag(C*n(2*nn+1:end))];
    ss = svd(A);
    singvals = [singvals, ss(end)];
    
    [U, V, W] = svd(A);
    q = W(:, end);
    if abs(sum(sign(q) - sign(q0))) > length(q)/2
        q = -q;
    end

    D = C'*diag(q)*C;
    [phi, lambda] = eig(D);
    B = phi(:, 1:4);

    n(1:nn) = B*pinv(B)*n(1:nn);
    n(nn+1:2*nn) = B*pinv(B)*n(nn+1:2*nn);
    n(2*nn+1:end) = B*pinv(B)*n(2*nn+1:end);
    
end

norm(n-n0)
plot(singvals);

xlens = C*n0(1:nn);
ylens = C*n0(nn+1:2*nn);
zlens = C*n0(2*nn+1:end);
lens0 = sum([xlens ylens zlens].^2, 2).^0.5;
xlens = C*n(1:nn);
ylens = C*n(nn+1:2*nn);
zlens = C*n(2*nn+1:end);
lens = sum([xlens ylens zlens].^2, 2).^0.5;
lengthchange = norm(lens0-lens)


% Display tensegrity:
nodes = [n(1:nn) n(nn+1:2*nn) n(2*nn+1:end)];
nodes0 = [n0(1:nn) n0(nn+1:2*nn) n0(2*nn+1:end)];
nodes(:,3) = nodes(:,3) - min(nodes(:,3)); % Make minimum node z=0 height.
nodes0(:,3) = nodes0(:,3) - min(nodes0(:,3)); % Make minimum node z=0 height.

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
superBallplot2 = TensegrityPlot(nodes0, strings, bars, barRad*2, stringRad);
f = figure('units','normalized', 'outerposition',[0 0 1 1]);
generatePlot(superBallplot, gca);
generatePlot(superBallplot2, gca);
updatePlot(superBallplot);
updatePlot(superBallplot2);
%settings to make it pretty
axis equal
view(90, 0); % X-Z view
%view(3)
grid on
light('Position',[0 0 10],'Style','local')
lighting flat
colormap([0.8 0.8 1; 0 1 1])
lims = l;
xlim([-lims lims])
ylim([-lims lims])
zlim(1.6*[-0.01 lims])

%% Check stability
A = [transpose(C) * diag(C*n(1:nn));
     transpose(C) * diag(C*n(nn+1:2*nn));
     transpose(C) * diag(C*n(2*nn+1:end))];
Ainv = pinv(A, 1e-3);
AinvA = Ainv * A;

% Full V.
V = (eye(length(AinvA)) - AinvA);
% Find an orthonormal basis for the range of V.
% NOTE: We don't use orth here, because orth uses precision
% eps, which is too small!
[Q,R,E] = qr(V);
[m , ~] = size(R);
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
q = q0;
if size(V,2) > 0 % Can only do this if there is nullspace!
    p = zeros(12*3,1);
    f = 2*transpose(p)*transpose(Ainv)*V;
    c = 10;
    options = optimoptions('quadprog','Algorithm',  'interior-point-convex','Display','final');
    w = quadprog(2*transpose(V)*V, f, -V(1:s, :), zeros(s,1) - c + Ainv(1:s,:)*p, ...
        [],[], [],[], [], options);

    % Now w has been found, we can find q:
    q = Ainv * p + V * w;
end

D = transpose(C)*diag(q)*C;
fprintf('Nullity(D) = %d\n', 12 - rank(D, 1e-3));

% Check if stable:
xlens = C*n(1:nn);
ylens = C*n(nn+1:2*nn);
zlens = C*n(2*nn+1:end);
actuallens = sum([xlens ylens zlens].^2, 2).^0.5;
K_string = 1000;
K_bar = 100e3;
Q = diag(q);
S = diag([K_string*ones(s,1); K_bar*ones(b,1)]);
Ke = A*diag(ones(length(actuallens),1)./actuallens)^2*(S-Q)*A';
Kg = kron(eye(3), C'*Q*C);
K = Ke + Kg;
eigs = sort(eig(K));
fprintf('Stiffness matrix eigenvalues:\n');
disp(eigs(1:10));
if norm(eigs(1:6)) < 1e-5 && all(real(eigs)>-1e-3)
    fprintf('Stable!\n');
else
    fprintf('Unstable!\n');
end