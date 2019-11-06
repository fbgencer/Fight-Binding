% create a 2D grid
th = linspace(0,pi,200);    % inclination
phi = linspace(0,2*pi,200); % azimuth
[th,phi] = meshgrid(th,phi);
% compute spherical harmonic of degree 3 and order 1
Y = harmonicY(2,-2,th,phi,'type','real');
% plot the magnitude
r = abs(Y);
[x,y,z] = sph2cart(phi,pi/2-th,r);

s = surf(x+10,y+10,z+10,r,'EdgeColor','None');
light               % create a light
lighting gouraud    % preferred method for lighting curved surfaces