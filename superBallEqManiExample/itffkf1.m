% Iterative form-find under known external load forces.
% Alex Popescu

clear all;
clc;

% Logging
feaserr = [];
options = optimoptions('quadprog','Algorithm',  'interior-point-convex','Display','final');

% External load forces:
nodalmass = 3; % kg
p = [zeros(24,1); [3 -1 3 -1 3 -1  -1 -1 -1 -1 -1 -1]'] * nodalmass * 9.81;

%% 0. Start with initial node positions n and valid force densities:
alpha = 60/180*pi;
delta = 55/180*pi;
delta = acos(1/sqrt(3));

l = 1.65;
b = sqrt(3/8)*l;
n0 = SVDB2nodes(alpha, delta, b, l);
n0 = n0(:); % Make column vector

% Connectivity matrix:
C = [0,0,0,0,0,1,0,0,-1,0,0,0;0,0,0,1,0,0,0,0,-1,0,0,0;0,0,0,1,0,0,-1,0,0,0,0,0;0,1,0,0,0,0,-1,0,0,0,0,0;0,1,0,0,0,0,0,0,0,0,-1,0;0,0,0,0,0,1,0,0,0,0,-1,0;1,0,0,0,0,-1,0,0,0,0,0,0;0,0,0,1,-1,0,0,0,0,0,0,0;0,1,-1,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,1,0,0,-1,0;0,0,0,0,0,0,1,0,0,-1,0,0;0,0,0,0,0,0,0,0,1,0,0,-1;1,0,0,0,0,0,0,0,0,0,-1,0;0,0,1,0,0,0,-1,0,0,0,0,0;0,0,0,0,1,0,0,0,-1,0,0,0;0,1,0,0,0,0,0,-1,0,0,0,0;0,0,0,1,0,0,0,0,0,-1,0,0;0,0,0,0,0,1,0,0,0,0,0,-1;1,0,0,0,-1,0,0,0,0,0,0,0;0,0,1,0,-1,0,0,0,0,0,0,0;1,0,-1,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,1,0,0,0,-1;0,0,0,0,0,0,0,0,0,1,0,-1;0,0,0,0,0,0,0,1,0,-1,0,0;1,-1,0,0,0,0,0,0,0,0,0,0;0,0,1,-1,0,0,0,0,0,0,0,0;0,0,0,0,1,-1,0,0,0,0,0,0;0,0,0,0,0,0,1,-1,0,0,0,0;0,0,0,0,0,0,0,0,1,-1,0,0;0,0,0,0,0,0,0,0,0,0,1,-1];
s = 24;
b = 6;
nn = 12; % Number of nodes

% Compute A ,the equilibrium matrix
A = [transpose(C) * diag(C*n0(1:nn));
     transpose(C) * diag(C*n0(nn+1:2*nn));
     transpose(C) * diag(C*n0(2*nn+1:end))];
n = n0;

% Ainv = pinv(A, 1e-3);
% AinvA = Ainv * A;
% 
% % Full V.
% V = (eye(length(AinvA)) - AinvA);
% % Find an orthonormal basis for the range of V.
% % NOTE: We don't use orth here, because orth uses precision
% % eps, which is too small!
% [Q,R,E] = qr(V);
% [m , n] = size(R);
% j=1;
% i=1;
% while i<=m
%     if norm(R(i,:))>10^-12
%         R_new(j,:)=R(i,:);
%         j=j+1;
%     else
%         i=m;
%     end
%     i=i+1;
% end
% V=Q(:,1:j-1);
% p = zeros(12*3,1);
% f = 2*transpose(p)*transpose(Ainv)*V;
% c = 10;
% options = optimoptions('quadprog','Algorithm',  'interior-point-convex','Display','off');
% w = quadprog(2*transpose(V)*V, f, -V(1:s, :), zeros(s,1) - c + Ainv(1:s,:)*p, ...
%     [],[], [],[], [], options);
% 
% % Now w has been found, we can find q_s:
% q = Ainv * p + V * w;
% 
% return;

Niter = 100;
for niter = 1:Niter
    %error_before = norm(A*q-p);
    %% 1. Compute feasible force densities
    %     q = pinv(A)*p;
    %     if any(q(1:s)<0)
    %          fprintf('q negative!\n');
    %          break;
    %     end
    xlens = C*n(1:nn);
    ylens = C*n(nn+1:2*nn);
    zlens = C*n(2*nn+1:end);
    lengths = sum([xlens ylens zlens].^2, 2).^0.5;
    qmin = ones(length(lengths),1)*10 ./ lengths; % min/max q limits for strings
    qmax = ones(length(lengths),1)*190 ./ lengths;
    qmin(end-b+1:end) = -1e6; % No limits on q in bars
    qmax(end-b+1:end) = 1e6;
    q = quadprog(2*(A')*A + 0*eye(length(lengths))*norm(A)/1000,... % H
                 -2*p'*A,... % f
                 [],... % Inequality constraint matrix
                 [],... % Inequality constraint vector
                 [],... % Equality constraint matrix
                 [],... % Equality constraint vector
                 qmin,... % q lower bound
                 qmax,... % q upper bound
                 [],... % x0
                 options);
             
    %q = [ones(s,1); -1.63*ones(b,1)];
    %cost = q'*(2*(A')*A + eye(length(q))*norm(A)/100)*q + (-2*p'*A)*q
    %error_after = norm(A*q-p);
    
    %% 2. Compute delta n to increase feasibility
    % Compute constraint matrix F(n):
    Cb = C(end-b+1:end, :); % Bar connectivity matrix is last b rows of C.
    F = kron([diag(Cb*n(1:nn)) diag(Cb*n(nn+1:2*nn)) diag(Cb*n(2*nn+1:end))], [1, -1]);
    % Compute force density matrix D(q):
    D_1d = transpose(C)*diag(q)*C;
    eig(D_1d)
    [V,E] = eig(D_1d); % Eigenvalue decomposition, D_1d=V*E*V'
    E(1,1) = 0;
    E(2,2) = 0;
    E(3,3) = 0;
    E(4,4) = 0;
    D_1d = V*E*V';
    
    D = kron(eye(3), D_1d); % Make 3n x 3n D matrix.
    dneps = 0.01; % "small" value of delta n_i is around 1cm
    fixedNodes = repmat(logical([1 0 1 0 1 0 zeros(1,6)]), [1 3]);
    dnub = ones(length(n),1)*dneps;
    dnlb = -ones(length(n),1)*dneps;
    dnub(24+[1 3 5]) = 0; % Enforce delta z=0 for bottom nodes to prevent going up/down
    dnlb(24+[1 3 5]) = 0;
    tooHighNodes = [logical(zeros(24,1)); n(25:end)>max(n0(25:end))];
    dnub(tooHighNodes) = 0; % Prevent going up.
    %dnub(24+[8 10 12]) = 0; % Enforce delta z<0 for top nodes to prevent going up
    dn = quadprog(2*transpose(D)*D + eye(length(n))*norm(D)/100,... % H
                 2*(n'*(D')*D - p'*D),... % f
                 [],... % Inequality constraint matrix
                 [],... % Inequality constraint vector
                 F,... % Equality constraint matrix
                 zeros(b,1),... % Equality constraint vector
                 dnlb,... % dn lower bound
                 dnub,... % dn upper bound
                 [],... % x0
                 options);
    n = n + dn;
    
    % Compute the new A for new nodes:
    A = [transpose(C) * diag(C*n(1:nn));
         transpose(C) * diag(C*n(nn+1:2*nn));
         transpose(C) * diag(C*n(2*nn+1:end))];
     
     % SVD of A:
    [U,S,V] = svd(A);
    %S(30,30) = 0;
    A = U*S*V';
    sv0 = S(30,30);
    
    error = norm(A*q-p)
    feaserr(end+1) = sv0;
    
    
    %% 3. Enforce length of bars is exactly the same:
    xbl = Cb*n(1:nn);
    ybl = Cb*n(nn+1:2*nn);
    zbl = Cb*n(2*nn+1:end);        
    xbc = abs(Cb)*n(1:nn)/2;
    ybc = abs(Cb)*n(nn+1:2*nn)/2;
    zbc = abs(Cb)*n(2*nn+1:end)/2;
    for bb = 1:b
        n([bb*2-1; 12+bb*2-1; 24+bb*2-1]) = [xbc(bb); ybc(bb); zbc(bb)] + [xbl(bb); ybl(bb); zbl(bb)]/norm([xbl(bb); ybl(bb); zbl(bb)])*l/2;
        n([bb*2; 12+bb*2; 24+bb*2]) = [xbc(bb); ybc(bb); zbc(bb)] - [xbl(bb); ybl(bb); zbl(bb)]/norm([xbl(bb); ybl(bb); zbl(bb)])*l/2;
    end
    
    % Check length constraint error:
    xbl = Cb*n(1:nn);
    ybl = Cb*n(nn+1:2*nn);
    zbl = Cb*n(2*nn+1:end);
    barlenerr = sum([xbl ybl zbl].^2, 2).^0.5 - ones(b,1)*l;
    
    [U,S,V] = svd(A);
    S(30,30) = 0;
    A = U*S*V';
end

% Check if stable:
xlens = C*n(1:nn);
ylens = C*n(nn+1:2*nn);
zlens = C*n(2*nn+1:end);
actuallens = sum([xlens ylens zlens].^2, 2).^0.5;
K_string = 1e3;
K_bar = 1000e3;
Q = diag(q);
S = diag([K_string*ones(s,1); K_bar*ones(b,1)]);
Ke = A*diag(ones(length(actuallens),1)./actuallens)^2*(S-Q)*A';
Kg = kron(eye(3), C'*Q*C);
K = Ke + Kg;
eigs = eig(K);
eigs(end-10:end) % Show a few eigenvalues

plot(feaserr);
%return;

% Calculate rest lengths:
actuallens = actuallens(1:s);
T_cmd = q(1:s) .* actuallens(1:s);
K = 1e3;
restlens = actuallens - T_cmd / K;

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
superBallplot2 = TensegrityPlot(nodes0, strings, bars, barRad, stringRad);
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