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

close all; clear; clc;
gp = lattice_drawer(figure(1),10,10,10);
gp.xaxis_symmetric();
gp.yaxis_symmetric();
gp.zaxis_symmetric();
grid;


a = 3;
% gp.draw('line',0,0,0,a,0,0);
% gp.draw('line',a,0,0,a,a,0);
% gp.draw('line',a,a,0,0,a,0);
% gp.draw('line',0,a,0,0,0,0);
% gp.draw('line',0,0,0,0,0,a);
% gp.draw('line',0,0,a,a,0,a);
% gp.draw('line',a,0,a,a,a,a);
% gp.draw('line',a,a,a,0,a,a);
% gp.draw('line',0,a,a,0,0,a);
% gp.draw('line',0,a,a,0,0,a);
% gp.draw('line',a,0,0,a,0,a);
% gp.draw('line',a,a,0,a,a,a);
% gp.draw('line',0,a,0,0,a,a);




cc = gp.draw('cuboid',0,0,0,a,a,a);
%cc = gp.draw('cuboid',a,a,a,a,a,a);
%sp = gp.draw('sphere',0,0,0,0.2);
atoma = gp.draw('sphere blue',0,0,0,0.25,'Visible','off');
atomc = gp.copy_to(atoma,0,a,a,'Visible','on');
%hold on;
%[x,y,z] = sphere();
%fvc = surf2patch(2*x-2,2*y-3,2*z,z); 
%p = patch(gp.fig.CurrentAxes,fvc,'LineStyle','None','FaceColor','blue'); 
%shading faceted; 






%o1 = gp.draw('line',3,3,3,10);
%o1 = gp.draw('line',3,10,25,20);
%o1 = gp.draw('line',25,20,3,3);

% a = 3;
% a1 = [a*sqrt(3)/2, a/2,0];
% a2 = [a*sqrt(3)/2, -a/2,0];

% gp.set_title('$Graphene\hspace{1mm}C_{6}$','Interpreter','latex');
% gp.set_xlabel('$X$','Interpreter','latex');
% %x = 5; y = 5;

% cx = 15; cy = 15;

% for i =  0:1
% for j = 0:1
    
%     cx = i*a1(1)+j*a2(1)+cx;
%     cy = i*a1(2)+j*a2(2)+cy;
    
%     l1 = gp.draw('line black',cx-1,cy,cx+1,cy);
%     cr = gp.draw('circle red',cx-1,cy,0.5);
%     cb = gp.draw('circle blue',cx+1,cy,0.3,'Visible','off');
%     cx = 15; cy = 15;
% end
% end

% dumy = gp.draw('rect',5,5,3,4,'FaceColor','blue','Visible','off');
% a = gp.copy_to(dumy,10,25);
% a.Visible = 'on';

%o1 = gp.draw('circle black',5,5,3);
%o2 = gp.copy_to(o1,25,25,1);


%Working code
% gp.fig.Position = [713 175 560 420];
% ctr = 1;
% for i =  -7:7
% for j =  -7:7
%     
%     cx = i*a1(1)+j*a2(1)+cx;
%     cy = i*a1(2)+j*a2(2)+cy;
% 
% cxa= cx-1;
% cya = cy;
% cxb= cx+1;
% cyb = cy;
% 
% 
% %next first neighbor
% cx1 = cx+a1(1);
% cy1 = cy+a1(2);
% cx1a = cx1-1;
% cy1a = cy1;
% cx1b = cx1+1;
% cy1b = cy1;
% l1 = gp.draw('line',cxb,cyb,cx1a,cy1a);
% 
% cx1 = cx-a1(1);
% cy1 = cy-a1(2);
% cx1a = cx1-1;
% cy1a = cy1;
% cx1b = cx1+1;
% cy1b = cy1;
% l1 = gp.draw('line',cxa,cya,cx1b,cy1b);
% 
% cx1 = cx-a2(1);
% cy1 = cy-a2(2);
% cx1a = cx1-1;
% cy1a = cy1;
% cx1b = cx1+1;
% cy1b = cy1;
% l1 = gp.draw('line',cxa,cya,cx1b,cy1b);
% 
% cx1 = cx+a2(1);
% cy1 = cy+a2(2);
% cx1a = cx1-1;
% cy1a = cy1;
% cx1b = cx1+1;
% cy1b = cy1;
% l1 = gp.draw('line',cxb,cyb,cx1a,cy1a);
%     
%     l1 = gp.draw('line black',cx-1,cy,cx+1,cy);
%     cr = gp.draw('circle',cx-1,cy,0.3);
%     cb = gp.draw('circle blue',cx+1,cy,0.3);
%     %gp.set_text(cr,sprintf('%d',ctr));
%     cx = 15; cy = 15;
% end
% end







%gp.infinite_layer();

% for j = 0:3*a:20
% for i = 0:6
%     h1=gp.draw('hexagon',j+x,y+i*(sqrt(3)*a),a);
%     h2=gp.draw('hexagon',j+x+3*a/2,a*sqrt(3)/2+y+i*(sqrt(3)*a),a);
%     for kk = 1:6
%         ll = h1.UserData(kk);
%         ll2 = h2.UserData(kk);
%         ll.Marker = 'o';
%         ll.MarkerFaceColor = 'blue';
%         ll2.Marker = 'o';
%         ll2.MarkerFaceColor = 'blue';
%     end
%     gp.draw('circle',j+x,y+i*(sqrt(3)*a),0.3,'black');
%     gp.draw('circle',j+x+3*a/2,a*sqrt(3)/2+y+i*(sqrt(3)*a),0.3,'green');
% end
% end

%cc = gp.draw('circle',10,20,5);
%sv = gp.draw('line',10,0,10,20);


%gp.draw('line',35,0,35,40);
%gp.draw('line',0,35,40,35);
%v1=gp.draw('vector',30,30,35,35,'red');

%v2=gp.draw('vector',30,30,38,33,'black');
%v3=gp.draw('vector',0,0,30,40,'black');
%gp.set_text(v1,'a1');
%gp.set_text(v3,'XYZ');

%gp.save_pdf('hexa_lattice2');

%gp.draw('line',5,5,20,5);
% 4

% gp.draw('hexagon',10,10+5*sqrt(3),5);
% gp.draw('hexagon',10+15/2,10+5*sqrt(3)/2,5);
% gp.draw('hexagon',10+15/2,10-5*sqrt(3)/2,5);
% gp.draw('circle',10,10,1);
% gp.draw('circle',10,10+5*sqrt(3),1);
% gp.draw('circle',10+15/2,10+5*sqrt(3)/2,1);
% gp.draw('circle',10+15/2,10-5*sqrt(3)/2,1);
% 
% gp.draw('vector',10,10,10,10+5*sqrt(3),'black');
% gp.draw('vector',10,10,10+15/2,10+5*sqrt(3)/2,'black');
% gp.draw('vector',10,10,10+15/2,10-5*sqrt(3)/2,'black');



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


%1- Atomic species, positions and orbitals verince unit cell'i çıkartsın ve unit cell'in 2D veya 3D interaktif görüntüsünü çıkartabilsin
%2- bantları çıkarttıktan sonra electronic DOS çıkartılması (toplam ve her band için ayrı ayrı istenebilsin)
%3- Brillouin Zone hesaplasın, reciprocal space vektörlerini bulsun, Brillouin zone'daki yüksek simetrisi olan noktaları çıkartsın
%3- elektron bantların, orbital, angular simetri, spin ve atomic ağırlıklarını versin
%4- 2D veya 3D Fermi yüzeyini çıkartsın, içeri doğru topolojisini versin, Fermi wave-vectorleri versin
%5- introducing the spin and the spin-orbit coupling
%6- from tight binding to electron-phonon interaction
