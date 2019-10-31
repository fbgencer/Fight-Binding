close all; clear; clc;
Eo = 0;
t = 0.8;
a0 = 1.42*1e-10;
a = a0*sqrt(3);

from = 2*pi/a;
len = 300;
kx = linspace(-from,from,len);
[kx,ky] = meshgrid(kx);
E = Eo+t.*sqrt(1+4*cos( (sqrt(3)/2)*kx*a).*cos(ky*a/2)+4*((cos(ky*a/2).^2)));
%E = Eo + t.*sqrt(E);

f = figure(1);
surf(kx,ky,E);
hold on;
surf(kx,ky,-E);
f = figure(2);
contour(kx,ky,E);
%colormap(f,jet);

%Gamma = [0, 0]
%K1 = [-4*pi / (3*sqrt(3)*a_cc), 0]
%M = [0, 2*pi / (3*a_cc)]
%K2 = [2*pi / (3*sqrt(3)*a_cc), 2*pi / (3*a_cc)]


Gamma = [0, 0];
K1 = [-2*pi / (sqrt(3)*a), 2*pi / (3*a)];
M = [2*pi / (sqrt(3)*a), 0];
K2 = [2*pi / (sqrt(3)*a), 2*pi / (3*a)];

%kx=(2*pi/(sqrt(3)*a));
%ky=(2*pi/(3*a))*s ;

f2 = figure(4);
kx = [linspace(K1(1),Gamma(1),len),linspace(Gamma(1),M(1),len),linspace(M(1),K2(1),len)];
ky = [linspace(K1(2),Gamma(2),len),linspace(Gamma(2),M(2),len),linspace(M(2),K2(2),len)];

%E1=t.*sqrt(1+4.*cos((sqrt(3).*a/2).*kx).*cos((a/2).*ky)+4.*cos((a/2).*ky).^2) ;
%E2=-E1; 

%kx = [linspace(Gamma(1),M(1),len),linspace(M(1),K2(1),len),linspace(K2(1),Gamma(1),len)];
%ky = [linspace(Gamma(2),M(2),len),linspace(M(2),K2(2),len),linspace(K2(2),Gamma(2),len),];


Esym = Eo-t.*sqrt(1+4.*cos( (sqrt(3)/2).*kx.*a).*cos(ky.*a/2)+4.*((cos(ky.*a/2).^2)));
Esymp = Eo+t.*sqrt(1+4.*cos( (sqrt(3)/2).*kx.*a).*cos(ky.*a/2)+4.*((cos(ky.*a/2).^2)));
plot(Esym)
hold on;
plot(Esymp)
