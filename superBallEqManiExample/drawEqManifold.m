% Draw the Equilibrium Manifold parameterized by alpha, delta, b:

npts = 100;

% Fixed parameters:
l = 1;

%cum_space = zeros(npts, npts);

%for b = linspace(0, 1.5, 50)
for b = 0.2188/1.65
    space = zeros(npts, npts); % 2d space for (alpha, delta)
    
    alpha_i = linspace(0, 90, npts);
    delta_i = linspace(0, 90, npts);
    
    % Plot boundary lines:
%     c = b./(2*sqrt(3)*l*sin(delta_i/180*pi));
%     acosc = acos(c);
%     acosc(imag(acosc)~=0) = NaN;
    %plot(delta_i, (acosc-pi/6)/pi*180); hold on;
    %plot(delta_i, (-acosc-pi/6 + pi)/pi*180); hold on;
    %alpha2 = (acos(b./(2*sqrt(3)*l*sin(delta_i/180*pi))) - pi/6) * 180/pi;
    %alpha1(imag(alpha1)~=0) = NaN;
    %alpha2(imag(alpha2)~=0) = NaN;
    %plot(delta_i, alpha1, delta_i, alpha2);
    
    %return;
    [A,D] = meshgrid(alpha_i, delta_i);
    Z = A*0 + NaN;
    S = A*0 + NaN;
    
    
        for x = 1:npts
            for y = 1:npts
                alpha = alpha_i(x)/180*pi;
                delta = delta_i(y)/180*pi;
                if l*sin(delta)*abs(cos(alpha + pi/6)) < b/(2*sqrt(3)) ...
                && sin(alpha + pi/6) < 3*l*sin(delta)/(2*b)
                    Z(x, y) = b;
                    u = sin(delta) * cos(alpha + pi/6.0);
h = cos(delta)/(2.0*u) * (l*u + sqrt(b*b/3.0 - 3.0*l*l*u*u) - b/sqrt(3.0));
                    % Saddle lengths
                    S(x, y) = sqrt(h*h + b*b/3.0 + l*l*sin(delta)*sin(delta) ...
           - 2.0/sqrt(3.0)*l*b*sin(delta)*cos(alpha - pi/6.0) );
                end
            end
        end
    
    s = surf(A, D, Z);
    %hold on;
    %s.EdgeColor = 'none';
    %surf(A, D, S);
    %cum_space(boolean(space)) = b;
end

%imagesc([0, delta_i(end)/pi*180], [0, alpha_i(end)/pi*180], cum_space);

%set(gca,'YDir','normal');
ylabel('alpha'); % y axis = axis 1
xlabel('delta'); % x axis = axis 2
zlabel('b');