close all; clear; clc;
Eo = 1;
t = 2;
a0 = 1.42*1e-10;
a = a0*sqrt(3);

from = 4*pi/a;

kx = linspace(-from,from,100);
[KX,KY] = meshgrid(kx);
E = abs(1+4*cos( (sqrt(3)/2)*KX*a).*cos(KY*a/2)+4*(cos(KY*a/2).*cos(KY*a/2)));
E = Eo + t.*sqrt(E);

f = figure(1);
surf(KX,KY,E);
hold on;
surf(KX,KY,-E);
colormap(f,jet);

K1 = -4*pi / (3*sqrt(3)*a);
M = 2*pi / (3*a);
K2 = 2*pi / (3*sqrt(3)*a);
Gamma = 0;
xticks([K1 Gamma K2 M]);
xticklabels({'K','K','M','G'});
%,'M','K'
