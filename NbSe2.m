% NbSe2
clear all; clc; close all;
a = 1;

%primitive vectors
a1 = [3.491454 0.000000 0.000000];
a2 = [-1.745727 3.023688 0.000000];
a3 = [0.000000 0.000000 20.902414];

tb = tightbinding(3,a1,a2,a3);% Start with name and primitive vectors
% Nb1 = [0.666667 0.333333 0.000807];
% Nb2 = [0.333333 0.666667 0.334141];
% Nb3 = [0.000000 0.000000 0.667474];
% Se1 = [0.333333 0.666667 0.081375];
% Se2 = [0.000000 0.000000 0.253484];
% Se3 = [0.000000 0.000000 0.414709];
% Se4 = [0.666667 0.333333 0.586817];
% Se5 = [0.666667 0.333333 0.748042];
% Se6 = [0.333333 0.666667 0.920151];

Nb1 = [1.745729    1.007895    0.0168689];
Nb2 = [-0.000002    2.015793    6.9843549];
Nb3 = [ 0.000000    0.000000    13.951818];
Se1 = [ -0.000002    2.015793    1.700934];
Se2 = [ 0.000000    0.000000    5.298428];
Se3 = [ 0.000000    0.000000    8.668420];
Se4 = [ 1.745729    1.007895    12.265893];
Se5 = [ 1.745729    1.007895    15.635883];
Se6 = [ -0.000002    2.015793    19.233377];

tb.set_unit_cell('Nb',Nb1,'Nb',Nb2,'Nb',Nb3,...
	'Se',Se1,'Se',Se2,'Se',Se3,'Se',Se4,'Se',Se5,'Se',Se6);

tb.set_orbital('3s','px,py,pz');
tb.set_metric_unit('A');
tb.spin_orbit_coupling('false');
%tb.add_hrfile('NbSe2_hr.dat');
tb.set_fermi_level(8.6675);



if(0)
% 0.0000 0.0000 0.0000 20 !Gamma
% 0.5000 0.0000 0.0000 20 !M
% 0.3300 0.3300 0.0000 20 !K
% 0.0000 0.0000 0.0000 20 !Gamma
% 0.5000 0.0000 0.5000 20 !L
G = {[0 0 0],'$$\Gamma$$'};
M = {2*pi*[0.5 0 0],'M'};
K = {2*pi*[0.33 0.33 0],'K'};
L = {2*pi*[0.5 0 0.5],'L'};

precision = [100,30,100];

fig_hsym = figure("Name","High Symmetry Points Figure");
tb.plot_high_symmetry_points(fig_hsym,precision,G,M,K,G);
ax = fig_hsym.CurrentAxes;
ax.YLim = [-5 5];
ax.YTick = -5:1:5;
ax.YTickLabel = {'-5','-4','-3','-2','-1','0','1','2','3','4','5'}
%,W,L,K,G
end




if(1)
lat_f = figure("Name","Lattice Figure");
gp = lattice_drawer(lat_f,20,20,20);
gp.xaxis_symmetric();
gp.yaxis_symmetric();
gp.zaxis_symmetric();


%Gallium covalent radius = 1.26 ang, 1.36 atomic
%As covalent radius = 1.19 ang, 1.14 atomic

atoma = gp.draw('sphere rgb:4589e4',0,0,0,1.36/2,'Visible','off');
atomb = gp.draw('sphere black',0,0,0,1.14/2,'Visible','off');
bond = gp.draw('line red',0,0,0,0,'Visible','off','LineWidth',4,'LineStyle','-');
%rgb:FF55FF
bonds = {bond};
atoms = {atoma,atomb};
from_to = 10:40;
tb.plot_lattice(gp,"bonds",bonds,"atoms",atoms);



%an = a;
%gp.draw('cuboid',an,an,an,an,an,an,'FaceColor','None','LineWidth',1);

end