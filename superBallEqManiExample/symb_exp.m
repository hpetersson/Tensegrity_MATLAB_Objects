% Symbolic experiments
clear all;

% Connectivity matrix.
C = [0,0,0,0,0,1,0,0,-1,0,0,0;0,0,0,1,0,0,0,0,-1,0,0,0;0,0,0,1,0,0,-1,0,0,0,0,0;0,1,0,0,0,0,-1,0,0,0,0,0;0,1,0,0,0,0,0,0,0,0,-1,0;0,0,0,0,0,1,0,0,0,0,-1,0;1,0,0,0,0,-1,0,0,0,0,0,0;0,0,0,1,-1,0,0,0,0,0,0,0;0,1,-1,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,1,0,0,-1,0;0,0,0,0,0,0,1,0,0,-1,0,0;0,0,0,0,0,0,0,0,1,0,0,-1;1,0,0,0,0,0,0,0,0,0,-1,0;0,0,1,0,0,0,-1,0,0,0,0,0;0,0,0,0,1,0,0,0,-1,0,0,0;0,1,0,0,0,0,0,-1,0,0,0,0;0,0,0,1,0,0,0,0,0,-1,0,0;0,0,0,0,0,1,0,0,0,0,0,-1;1,0,0,0,-1,0,0,0,0,0,0,0;0,0,1,0,-1,0,0,0,0,0,0,0;1,0,-1,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,1,0,0,0,-1;0,0,0,0,0,0,0,0,0,1,0,-1;0,0,0,0,0,0,0,1,0,-1,0,0;1,-1,0,0,0,0,0,0,0,0,0,0;0,0,1,-1,0,0,0,0,0,0,0,0;0,0,0,0,1,-1,0,0,0,0,0,0;0,0,0,0,0,0,1,-1,0,0,0,0;0,0,0,0,0,0,0,0,1,-1,0,0;0,0,0,0,0,0,0,0,0,0,1,-1];

% Bar parametrization:
bx = sym('bx', [6, 1], 'real'); % bar base positions
by = sym('by', [6, 1], 'real');
bz = sym('bz', [6, 1], 'real');

alpha = sym('alpha', [6, 1], 'real');
delta = sym('delta', [6, 1], 'real');

L = 1;
% nx = [  bx(1);
%         bx(1) + L*sin(delta(1))*cos(alpha(1));
%         bx(2);
%         bx(2) + L*sin(delta(2))*cos(alpha(2));
%         bx(3);
%         bx(3) + L*sin(delta(3))*cos(alpha(3));
%         bx(4);
%         bx(4) + L*sin(delta(4))*cos(alpha(4));
%         bx(5);
%         bx(5) + L*sin(delta(5))*cos(alpha(5));
%         bx(6);
%         bx(6) + L*sin(delta(6))*cos(alpha(6)) ];
% ny = [  by(1);
%         by(1) + L*sin(delta(1))*sin(alpha(1));
%         by(2);
%         by(2) + L*sin(delta(2))*sin(alpha(2));
%         by(3);
%         by(3) + L*sin(delta(3))*sin(alpha(3));
%         by(4);
%         by(4) + L*sin(delta(4))*sin(alpha(4));
%         by(5);
%         by(5) + L*sin(delta(5))*sin(alpha(5));
%         by(6);
%         by(6) + L*sin(delta(6))*sin(alpha(6)) ];
% nz = [  bz(1);
%         bz(1) + L*cos(delta(1));
%         bz(2);
%         bz(2) + L*cos(delta(2));
%         bz(3);
%         bz(3) + L*cos(delta(3));
%         bz(4);
%         bz(4) + L*cos(delta(4));
%         bz(5);
%         bz(5) + L*cos(delta(5));
%         bz(6);
%         bz(6) + L*cos(delta(6)); ];
% The 36x1 node vector is now parametrized by 3*6+6*2=30 variables.

nx = sym('nx', [12,1], 'real');
ny = sym('ny', [12,1], 'real');
nz = sym('nz', [12,1], 'real');

x = sym('x', [30,1], 'real');

% Compute A, the equilibrium matrix, 36x30.
A = [transpose(C) * diag(C*nx);
     transpose(C) * diag(C*ny);
     transpose(C) * diag(C*nz)]
% 
% q = sym('