close all;
clear;
clc;
f = @(x,y,z) x.^2 + y.^2 + z.^2 - 5;
interval = [-5 5 -5 5 -2 5];
%fimplicit3(f,interval)

a = 1e-10;
from = -2*pi./a;
to = -from;
len = 100;
k = linspace(from,to,len);
[kx,ky,kz] = meshgrid(k);
Eo = 0;
t = 0.5;
En = -Eo - 4.*t.*(cos(ky.*a/2).*cos(kz.*a/2)+cos(kx.*a/2).*cos(kz.*a/2)+cos(kx.*a/2).*cos(ky.*a/2));
Ef1 = -Eo - 8*t;

E = @(kx,ky,kz) -Eo - 4.*t.*(cos(ky.*a/2).*cos(kz.*a/2)+cos(kx.*a/2).*cos(kz.*a/2)+cos(kx.*a/2).*cos(ky.*a/2)) + Eo - 1*t;
int = [from to from to from to];
fimplicit3(E,int);

%contour3(kx(:,:,1),ky(:,:,1),E(:,:,1))
%surf(kx(:,:,1),ky(:,:,1),E(:,:,1))

[kx,ky] = meshgrid(k);
Esq = -Eo - 2.*t.*(cos(kx.*a)+cos(ky.*a));
figure();
contour(kx,ky,Esq);
