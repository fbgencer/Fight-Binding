% NbSe2
clear all; clc; close all;
a = 1;

%primitive vectors
a1 = [3.493174    0.000000    0.000000];
a2 = [-1.746587    3.025177    0.000000];
a3 = [0.000000    0.000000    28.215334];

tb = tightbinding(3,a1,a2,a3);% Start with name and primitive vectors


Nb1 =[    0.000000    0.000000    0.000000];
Nb2 =[    0.000000    0.000000    7.062044];
Nb3 =[    -0.000002    2.016786    14.107667];
Nb4 =[    0.000000    0.000000    21.153290];
Se1 =[    -0.000002    2.016786    26.531019];
Se2 =[    1.746589    1.008391    19.465731];
Se3 =[    1.746589    1.008391    8.749603];
Se4 =[    -0.000002    2.016786    1.684314];
Se5 =[    1.746589    1.008391    22.837660];
Se6 =[    1.746589    1.008391    5.377674];
Se7 =[    0.000000    0.000000    12.421914];
Se8 =[    0.000000    0.000000    15.793421];


tb.set_unit_cell('Nb',Nb1,'Nb',Nb2,'Nb',Nb3,'Nb',Nb4,...
	'Se',Se1,'Se',Se2,'Se',Se3,'Se',Se4,'Se',Se5,'Se',Se6,'Se',Se7,'Se',Se8);

tb.set_orbital('px,py,pz','d1,d2,d3,d4,d5');
tb.set_metric_unit('m');
tb.spin_orbit_coupling('false');
tb.add_hrfile('NbSe2_hr_new.dat');
tb.set_fermi_level(8.5136);



if(1)
% 0.0000 0.0000 0.0000 20 !Gamma
% 0.5000 0.0000 0.0000 20 !M
% 0.3300 0.3300 0.0000 20 !K
% 0.0000 0.0000 0.0000 20 !Gamma
% 0.5000 0.0000 0.5000 20 !L
G = {[0 0 0],'$$\Gamma$$'};
M = {2*pi*[0.5 0 0],'M'};
K = {2*pi*[0.33 0.33 0],'K'};
L = {2*pi*[0.5 0 0.5],'L'};

precision = [100,60,100];

fig_hsym = figure("Name","High Symmetry Points Figure");
tb.plot_high_symmetry_points(fig_hsym,precision,G,M,K,G);
ax = fig_hsym.CurrentAxes;
%ax.YLim = [-5 5];
%ax.YTick = -5:1:5;
%ax.YTickLabel = {'-5','-4','-3','-2','-1','0','1','2','3','4','5'};
%,W,L,K,G
end




if(0)
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