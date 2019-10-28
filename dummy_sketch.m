close all
clear all
clc
t=-2.550; 
acc=1.44e-10; 
a=1.732*acc;

s=1:-0.02:0 ;
%kx=(2*pi/(sqrt(3)*a))*s ;
kx=linspace(2*pi/(sqrt(3)*a),0,size(s,2));
ky=linspace(2*pi/(3*a),0,size(s,2));
E1=t.*sqrt(1+4.*cos((sqrt(3).*a/2).*kx).*cos((a/2).*ky)+4.*cos((a/2).*ky).^2) ;
E2=-E1; 


plot(E1,'r') 
hold on 
plot(E2,'r')

E = [E1];
En = E2;

s=0:0.02:1;
kx=(2*pi/(sqrt(3)*a))*s;
ky=0;
E1=t.*sqrt(1+4.*cos((sqrt(3).*a/2).*kx).*cos((a/2).*ky)+4.*cos((a/2).*ky).^2);
E2=-E1;

%plot(E1,'*')
%hold on
%plot(E2,'*')

E = [E E1];
En = [En E2];

s=0:0.02:1;
kx=(2*pi/(sqrt(3)*a));
ky=(2*pi/(3*a))*s ;
E1=t.*sqrt(1+4.*cos((sqrt(3).*a/2).*kx).*cos((a/2).*ky)+4.*cos((a/2).*ky).^2);
E2=-E1; 

%plot(E1,'*')
%hold on
%plot(E2,'*')


E = [E E1];
En = [En E2];

plot(E)
hold on;
plot(En);


