function restLengths = eqManiSVDBRestLengths( alpha, delta, l, b, k_stiffness )
% Generate rest lengths according to equilibrium manifold open-loop control
% law.
    
% Compute the control SVDB tensions from Skelton, 2001:
u = sin(delta) * cos(alpha + pi/6.0);
h = cos(delta)/(2.0*u) * (l*u + sqrt(b*b/3.0 - 3.0*l*l*u*u) - b/sqrt(3.0));

S = sqrt(h*h + b*b/3.0 + l*l*sin(delta)*sin(delta) ...
           - 2.0/sqrt(3.0)*l*b*sin(delta)*cos(alpha - pi/6.0) );
V = sqrt(b*b + l*l - 2.0*l*b*sin(delta)*sin(alpha + pi/6.0));
D = sqrt(h*h + b*b/3.0 + l*l - 2.0/sqrt(3.0)*l*b*sin(delta)*sin(alpha) - 2.0*l*h*cos(delta));
B = b;

T_V = V/D * 1.0 / (sqrt(3.0)*cos(alpha+pi/6.0))*((l*cos(delta)/h - 1.0)*sin(alpha-pi/6.0) - cos(alpha));
T_S = S/D*(l*cos(delta)/h - 1.0);
T_B = 1.0 / (6.0*D) * (3.0*l*l*sin(delta)*cos(delta) + 6*b*h*cos(alpha-pi/3.0) - 6*l*h*sin(delta) - 2.0*sqrt(3.0)*b*l*cos(delta)*sin(alpha)) / (sqrt(3.0)*h*cos(alpha+pi/6.0));
T_D = 1.0;


% Handle discontinuity at alpha=pi/3
if abs(alpha - pi/3.0) < 0.001
    h = l * cos(delta)/2.0;
    S = sqrt(h*h + b*b/3.0 + l*l*sin(delta)*sin(delta) ...
           - 2.0/sqrt(3.0)*l*b*sin(delta)*cos(alpha - pi/6.0) );
    V = sqrt(b*b + l*l - 2.0*l*b*sin(delta)*sin(alpha + pi/6.0));
    D = sqrt(h*h + b*b/3.0 + l*l - 2.0/sqrt(3.0)*l*b*sin(delta)*sin(alpha) - 2.0*l*h*cos(delta));
    B;
    
    T_V = V/D * (3*l/(2*b)*sin(delta) - 1.0);
    T_S = 1;
    T_B = 1.0/(6.0*D) * (2.0*b*b - 9.0*l*b*sin(delta) + 9.0*l*l*sin(delta)*sin(delta))/b;
    T_D = 1;
end

% Tension scaling factor:
T_scaling = 1500/2 / norm([T_S, T_V, T_D, T_B]) / sqrt(6);
T_scaling = 230;
T_S = T_S * T_scaling
T_V = T_V * T_scaling
T_D = T_D * T_scaling
T_B = T_B * T_scaling

restLengths = [ones(6, 1) * (S - T_S/k_stiffness);
               ones(6, 1) * (V - T_V/k_stiffness);
               ones(6, 1) * (D - T_D/k_stiffness);
               ones(6, 1) * (B - T_B/k_stiffness) ]
%restLengths(restLengths < 0.1) = 0.1;

end

