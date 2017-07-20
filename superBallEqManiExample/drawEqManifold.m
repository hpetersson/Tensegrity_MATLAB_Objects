% Draw the Equilibrium Manifold parameterized by alpha, delta, b:

npts = 200;

% Fixed parameters:
l = 1;

%cum_space = zeros(npts, npts);

%for b = linspace(0, 1.5, 50)
for b = 0.5
    space = zeros(npts, npts); % 2d space for (alpha, delta)
    
    alpha_i = linspace(0, 90, npts);
    delta_i = linspace(0, 90, npts);
    [A,D] = meshgrid(alpha_i, delta_i);
    Z = A*0 + NaN;
    
        for x = 1:npts
            for y = 1:npts
                alpha = alpha_i(x)/180*pi;
                delta = delta_i(y)/180*pi;
                if l*sin(delta)*abs(cos(alpha + pi/6)) < b/(2*sqrt(3)) ...
                && sin(alpha + pi/6) < 3*l*sin(delta)/(2*b)
                    Z(x, y) = b;
                end
            end
        end
    
    s = surf(A, D, Z);
    hold on;
    s.EdgeColor = 'none';
    %cum_space(boolean(space)) = b;
end

%imagesc([0, delta_i(end)/pi*180], [0, alpha_i(end)/pi*180], cum_space);

%set(gca,'YDir','normal');
ylabel('alpha'); % y axis = axis 1
xlabel('delta'); % x axis = axis 2
zlabel('b');