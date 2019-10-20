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

% close all;
% xc = 0;
% yc = 0;
% a = 1;
% a3 = a*sqrt(3);
% pgon = polyshape(xc+[-a/2 a/2 a a/2 -a/2 -a],yc+[a3 a3 0 -a3 -a3 0]);
% xc = 2*a;
% yc = 0;
% a = 1;
% a3 = a*sqrt(3);
% pgon2 = polyshape(xc+[-a/2 a/2 a a/2 -a/2 -a],yc+[a3 a3 0 -a3 -a3 0]);
% 
% xc = a;
% yc = 3;
% a = 1;
% a3 = a*sqrt(3);
% pgon3 = polyshape(xc+[-a/2 a/2 a a/2 -a/2 -a],yc+[a3 a3 0 -a3 -a3 0]);
% 
% hold on;
% plot(pgon);
% plot(pgon2);
% plot(pgon3);
%%
%a = 3;
%pgon1 = nsidedpoly(6,'Center',[0 0],'SideLength',a);
%pgon2 = nsidedpoly(6,'Center',[a+a/2 a*sqrt(3)/2],'SideLength',a);
%pgon3 = nsidedpoly(6,'Center',[a+a/2 -a*sqrt(3)/2],'SideLength',a);
%pgon3 = rotate(pgon1,90);
%plot([pgon1 pgon2 pgon3]);

%axis([0 100 0 100])
%pos = [50 50 10 10]; 
%r = rectangle('Position',pos,'Curvature',[1 1]);
%r.FaceColor = 'red';
%r1 = rectangle('Position',[0 0 10 10]);

close all;
clear;
clc;
gp = lattice_drawer(1,30,30);

%o1 = gp.draw('line',3,3,3,10);
%o1 = gp.draw('line',3,10,25,20);
%o1 = gp.draw('line',25,20,3,3);

o = [];
for i = 0:5
    %o2 = gp.draw('eqtri',5+i*3,5,3);
    %o3 = gp.draw('circle',5+i*3,5,.1,'red');
    
    %gp.draw('hexagon',5+i*2,10,1);
    %gp.draw('circle',5+i*2,10,.1,'red');
    ;
end

%gp.draw('line',5,5,20,5);

gp.draw('hexagon',10,10,5);
gp.draw('hexagon',10,10+5*sqrt(3),5);
gp.draw('hexagon',10+15/2,10+5*sqrt(3)/2,5);
gp.draw('hexagon',10+15/2,10-5*sqrt(3)/2,5);
gp.draw('circle',10,10,1);
gp.draw('circle',10,10+5*sqrt(3),1);
gp.draw('circle',10+15/2,10+5*sqrt(3)/2,1);
gp.draw('circle',10+15/2,10-5*sqrt(3)/2,1);

gp.draw('vector',10,10,10+15/2,10+5*sqrt(3)/2,'black');
gp.draw('vector',10,10,10+15/2,10+5*sqrt(3)/2,'black');
gp.draw('vector',10,10,10+15/2,10+5*sqrt(3)/2,'black');

% o1 = gp.draw('circle',5,5,1,'green','green');
% o2 = gp.draw('circle',5,15,1,'red','red');
% o3 = gp.draw('circle',15,15,1,'blue','blue');
% o4 = gp.draw('circle',15,5,1,'black','black');
% 
% 
% gp.set_text(o1,'A','Color','red');
% gp.set_text(o2,'B','Color',[1 1 1]);
% gp.set_text(o3,'C','Color',[1 1 1]);
% gp.set_text(o4,'D','Color',[1 1 1]);
% %l.draw_rect(30,10,5,5,'red');
% gp.connect(o1,o2);
% gp.connect(o2,o3);
% gp.connect(o3,o4);
% gp.connect(o2,o4);
% gp.connect(o1,o4);


