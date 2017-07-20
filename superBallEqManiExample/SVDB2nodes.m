function nodes = SVDB2nodes( alpha, delta, b, l )
% Returns the node array for 2-stage SVDB structure.

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

end

