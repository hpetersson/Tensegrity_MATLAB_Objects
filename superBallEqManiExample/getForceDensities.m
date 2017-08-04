function q = getForceDensities( nodes, C, s, p )
% Returns valid force densities for all members (strings + bars), given
% some node positions nodes and external forces p.
% C = connectivity matrix
% s = number of strings

% Compute A ,the equilibrium matrix
nn = length(nodes)/3;
A = [transpose(C) * diag(C*nodes(1:nn));
     transpose(C) * diag(C*nodes(nn+1:2*nn));
     transpose(C) * diag(C*nodes(2*nn+1:end))];
 
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

% Minimum force density for strings.
c = 10;

f = 2*transpose(p)*transpose(Ainv)*V;
options = optimoptions('quadprog','Algorithm',  'interior-point-convex','Display','final');
w = quadprog(2*transpose(V)*V, f, -V(1:s, :), zeros(s,1) - c + Ainv(1:s,:)*p, ...
    [],[], [],[], [], options);
q = Ainv*p + V*w;

end