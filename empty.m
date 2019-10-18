% close all;clc;clear;
% a = 1.42*1e-10;
% a = 1;
% from = 4*pi/a;
% kx = linspace(-from,from,80);
% t = 5;
% figure(1);
% %x = -20:1:20;
% [KX,KY] = meshgrid(kx);
% e1 = 1;
% e2 = 2;
% ka = linspace(0,3,100);
% E = (e1+e2)+sqrt((e1-e2)^2+(8*t*t)*(1+cos(KX*a)));
% E = E/2;
% surf(KX,KY,E) 
% %plot(ka,E,ka,-E);
% %hold on
% %surf(KX,KY,-E)

close all;
xc = 0;
yc = 0;
a = 1;
a3 = a*sqrt(3);
pgon = polyshape(xc+[-a/2 a/2 a a/2 -a/2 -a],yc+[a3 a3 0 -a3 -a3 0]);
xc = 2*a;
yc = 0;
a = 1;
a3 = a*sqrt(3);
pgon2 = polyshape(xc+[-a/2 a/2 a a/2 -a/2 -a],yc+[a3 a3 0 -a3 -a3 0]);

xc = a;
yc = 3;
a = 1;
a3 = a*sqrt(3);
pgon3 = polyshape(xc+[-a/2 a/2 a a/2 -a/2 -a],yc+[a3 a3 0 -a3 -a3 0]);

hold on;
plot(pgon);
plot(pgon2);
plot(pgon3);