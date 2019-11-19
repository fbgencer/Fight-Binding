%Graphene Structure
% [A]--[B]         [A]--[B]
% (1)    \        /  (2)
%         [A]--[B] >>> this is 0
% (3)    /        \  (4)
% [A]--[B]         [A]--[B]

%From (0) to go (1) we have -a1, so we will add hopping as following
%Our translation vector R = -a1 + 0a2 so it is [-1 0]
%Hopping amplitude as t
%interaction is which atoms interact which ones, in our case 0 to 1 A atoms
%interacts with B
%add_hopping(amplitude,interaction between pairs,translation vec)
%add_hopping(t,1,2,[0 0]);  % inside 0 A and B interaction interaction
%between zero and one which A and B respectively.
%add_hopping(t,1,2,[-1 0]);  % 0 -> 1
%add_hopping(t,2,1,[1 0]);   % 0 -> 2
%add_hopping(t,2,1,[0 1]);  % 0 -> 3
%add_hopping(t,1,2,[0 -1]);  % 0 -> 4

clear; clc; close all;
%graphene consts

a0 = 1.42*1e-10;
a = a0*sqrt(3);

%primitive vectors
a1 = [a*sqrt(3)/2,  a/2,  0];
a2 = [a*sqrt(3)/2, -a/2, 0];

%a1 = [a,0,0];
%a2 = [a/2,a/2*sqrt(3),0];


tb = tightbinding(2,a1,a2);% Start with name and primitive vectors
tb.set_unit_cell('A',[-a0/2 0],'B',[a0/2 0]); %give unit cell atoms and their locations
%tb.set_unit_cell('A',[0, -a0/2],'B',[0, a0/2]); %give unit cell atoms and their locations


Eo = 0;
t = 0.1;
tb.add_hopping(Eo,1,1,[0 0]);  % 0 to 0
tb.add_hopping(Eo,2,2,[0 0]);  % 0 to 0

tb.add_hopping(-t,1,2,[0 0]);  % 0 to 0
%tb.add_hopping(-t,2,1,[0 0]);  % 0 to 0 % No need to include this anymore!

tb.add_hopping(-t,1,2,[-1 0]);  % 0 to 1 with -a1 vector
tb.add_hopping(-t,2,1,[1 0]);   % 0 to 2 with a1 vector
tb.add_hopping(-t,1,2,[0 -1]);  % 0 to 3
%tb.add_hopping(-t,2,1,[0 1]);  % 0 to 4

%Do this, just make sure to get symmetric hamiltonian
%tb.hermitian_hamiltonian();
if(1)
lat_f = figure("Name","Lattice Figure");
gp = lattice_drawer(lat_f,20,20);
atoma = gp.draw('circle rgb:0066cc',0,0,0.3,'Visible','off');
atomb = gp.draw('circle rgb:fb7100',0,0,0.3,'Visible','off');
bond = gp.draw('line black',0,0,0,0,'Visible','off','LineWidth',2);
bonds = {bond};
atoms = {atoma,atomb};
tb.plot_lattice(gp,"bonds",bonds,"atoms",atoms,"x",-5:5);

uc_rect{1} = gp.draw('line',-2,0,0,1,'Color','red','LineWidth',2);
uc_rect{2} = gp.draw('line',0,1,2,0,'Color','red','LineWidth',2);
uc_rect{3} = gp.draw('line',2,0,0,-1,'Color','red','LineWidth',2);
uc_rect{4} = gp.draw('line',0,-1,-2,0,'Color','red','LineWidth',2);
end


if(0)

range = 2*pi/a; 
precision = 100;
k = tb.set_kvector(-range,range,precision);
fig_band = figure("Name","Energy Band Figure");
surfaces = tb.plot_energy_band(fig_band,k,'surface','EdgeColor','None');

rp = lattice_drawer(fig_band);
lin = rp.draw('line rgb:FF8000',0,0,0,0,'Visible','off','ZData',0.5,'LineWidth',2);
tb.plot_brillouin_zone(rp,'plot points','plot lines');

%surfaces{1}.Visible = 'off';

view(0,90);


Gamma = [0, 0,0];
K1 = [-2*pi / (sqrt(3)*a), 2*pi / (3*a), 0];
M = [2*pi / (sqrt(3)*a), 0, 0];
K2 = [2*pi / (sqrt(3)*a), 2*pi / (3*a), 0];

fig_hsym = figure("Name","High Symmetry Points Figure");
tb.plot_high_symmetry_points(fig_hsym,K1,Gamma,M,K2);
title(fig_hsym.CurrentAxes,'High Symmetry Points','Interpreter','Latex');
xlabel(fig_hsym.CurrentAxes,'$$\Gamma K M K$$','Interpreter','Latex');
xticks(fig_hsym.CurrentAxes,0:precision:precision*4);
xticklabels(fig_hsym.CurrentAxes,{'K','\Gamma','M','K'});


band_dw = lattice_drawer(fig_band); 
band_dw.draw('vector red',K1(1),K1(2),0.5,Gamma(1),Gamma(2),0.5,'LineWidth',2,'MaxHeadSize',0.3);
band_dw.draw('vector black',Gamma(1),Gamma(2),0.5,M(1),M(2),0.5,'LineWidth',2,'MaxHeadSize',0.3);
band_dw.draw('vector green',M(1),M(2),0.5,K2(1),K2(2),0.5,'LineWidth',2,'MaxHeadSize',0.5);
l = band_dw.draw('vector black',K2(1),K2(2),0.5,K1(1),K1(2),0.5,'LineWidth',2,'MaxHeadSize',0.3);

t = band_dw.draw('text',K1(1),K1(2),0.6,'K','Interpreter','Latex','FontSize',13,'VerticalAlignment','bottom','HorizontalAlignment','right');
band_dw.draw('text',Gamma(1),Gamma(2),0.6,'$$\Gamma$$','Interpreter','Latex','FontSize',13,'VerticalAlignment','bottom','HorizontalAlignment','right');
band_dw.draw('text',M(1),M(2),0.6,'M','Interpreter','Latex','FontSize',13,'VerticalAlignment','bottom','HorizontalAlignment','right');
band_dw.draw('text',K2(1),K2(2),0.6,'K','Interpreter','Latex','FontSize',13,'VerticalAlignment','bottom','HorizontalAlignment','right');
end


%rp.draw('vector black',0,0,b1(1)+b2(1),b1(2)+b2(2));
%rp.draw('vector yellow',0,0,b1(1)-b2(1),b1(2)-b2(2));






%gp_byhand;
%f.Color = [0.5 0.5 0.2];

%K1 = [-4*pi / (3*sqrt(3)*a), 0,0];
%Gamma = [0, 0, 0];
%gp = lattice_drawer(f);
%ll = gp.draw('vector',K1(1),K1(2),Gamma(1),Gamma(2));
%ll.MaxHeadSize = 1;
%ll.LineWidth = 2;
%M = [0, 2*pi / (3*a),0];
%ll = gp.draw('vector',Gamma(1),Gamma(2),M(1),M(2));
%ll.MaxHeadSize = 1;
%ll.LineWidth = 2;
%K2 = [2*pi / (3*sqrt(3)*a), 2*pi / (3*a),0];
%ll = gp.draw('vector',M(1),M(2),K2(1),K2(2));
%ll.LineWidth = 2;
