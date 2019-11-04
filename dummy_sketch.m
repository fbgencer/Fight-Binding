% close all
% clear all
% clc
% t=-2.550; 
% acc=1.44e-10; 
% a=1.732*acc;
% 
% s=1:-0.02:0 ;
% %kx=(2*pi/(sqrt(3)*a))*s ;
% kx=linspace(2*pi/(sqrt(3)*a),0,size(s,2));
% ky=linspace(2*pi/(3*a),0,size(s,2));
% E1=t.*sqrt(1+4.*cos((sqrt(3).*a/2).*kx).*cos((a/2).*ky)+4.*cos((a/2).*ky).^2) ;
% E2=-E1; 
% 
% 
% plot(E1,'r') 
% hold on 
% plot(E2,'r')
% 
% E = [E1];
% En = E2;
% 
% s=0:0.02:1;
% kx=(2*pi/(sqrt(3)*a))*s;
% ky=0;
% E1=t.*sqrt(1+4.*cos((sqrt(3).*a/2).*kx).*cos((a/2).*ky)+4.*cos((a/2).*ky).^2);
% E2=-E1;
% 
% %plot(E1,'*')
% %hold on
% %plot(E2,'*')
% 
% E = [E E1];
% En = [En E2];
% 
% s=0:0.02:1;
% kx=(2*pi/(sqrt(3)*a));
% ky=(2*pi/(3*a))*s ;
% E1=t.*sqrt(1+4.*cos((sqrt(3).*a/2).*kx).*cos((a/2).*ky)+4.*cos((a/2).*ky).^2);
% E2=-E1; 
% 
% %plot(E1,'*')
% %hold on
% %plot(E2,'*')
% 
% 
% E = [E E1];
% En = [En E2];
% 
% plot(E)
% hold on;
% plot(En);


%3d plot genelleştirmesi için

close all;
clear; clc;
% s = scatter3(5,5,5);
% s.Marker  ='o';
% s.MarkerFaceColor = 'blue';
% hold on;
% s = scatter3(6,20,0);
% s.Marker  ='o';
% s.MarkerFaceColor = 'blue';
%hold on
%x = 0:0.01:2*pi;
%r = rectangle('Position',[1 2 5 6]);
%plot(x,sin(x));

%vert = [0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1];
%fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
%p = patch('Vertices',vert,'Faces',fac,'FaceVertexCData',hsv(6),'FaceColor','flat')
%,...'FaceVertexCData',hsv(6),'FaceColor','flat'

x1 = 0;
y1 = 2;
z1 = 0;
w = 5;
l = 10;
h = 1;


vertx = [x1, x1, x1, x1, x1+w, x1+w, x1+w, x1+w];
verty = [y1, y1+l, y1+l, y1, y1, y1+l, y1+l, y1];
vertz = [z1, z1, z1+h, z1+h, z1, z1, z1+h, z1+h];
fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
p = patch(vertx,verty,vertz,'red','Faces',fac);

line([x1, 5],[y1 5],[z1 5]);
%patch(x,y,z,'red');
