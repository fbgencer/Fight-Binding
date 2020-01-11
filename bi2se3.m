clear all;
clc; close all;

a1 = [-2.069  -3.583614  0.000000];
a2 = [2.069  -3.583614  0.000000];
a3 = [0.000   2.389075  9.546667];

%bu cinlinin yaptigi degil nette buldugum 541837'nin sabitleri

%a1 =   [4.191366    0.000000    0.000000];
%a2 =   [-2.095683    3.629830    0.000000];
%a3 =   [0.000000    0.000000    29.928688];

tb = tightbinding(3,a1,a2,a3);% Start with dimension and primitive vectors


%pos_bi1 = [0.3990    0.3990    0.6970];
%pos_bi2 = [0.6010    0.6010    0.3030];
%pos_se1 = [0.0000    0.0000    0.5000];
%pos_se2 = [0.2060    0.2060    0.1180];
%pos_se3 = [0.7940    0.7940    0.8820];

pos_bi1 = [0.000000    -1.194539    6.654027];
pos_bi2 = [0.000000    -3.583614    2.892640];
pos_se1 = [0.000000     1.194538    4.773334];
pos_se2 = [0.000000    -1.194538    1.126507];
pos_se3 = [0.000000    -3.583615    8.420160];



% posbi1=[-0.000002    2.419888    1.948388];
% posbi2=[-0.000002    2.419888    8.027831];
% posbi3=[0.000000    0.000000    11.924637];
% posbi4=[0.000000    0.000000    18.004051];
% posbi5=[2.095685    1.209942    21.900856];
% posbi6=[2.095685    1.209942    27.980299];
% posse1=[0.000000    0.000000    0.000000];
% posse2=[0.000000    0.000000    6.456067];
% posse3=[2.095685    1.209942    3.520153];
% posse4=[2.095685    1.209942    9.976219];
% posse5=[2.095685    1.209942    16.432287];
% posse6=[-0.000002    2.419888    13.496402];
% posse7=[-0.000002    2.419888    19.952469];
% posse8=[-0.000002    2.419888    26.408535];
% posse9=[0.000000    0.000000    23.472622];


%tb.set_unit_cell('Bi',posbi1,'Bi',posbi2,'Bi',posbi3,'Bi',posbi4,'Bi',posbi5,'Bi',posbi6,...
%	'Se',posse1,'Se',posse2,'Se',posse3,'Se',posse4,'Se',posse5,'Se',posse6,'Se',posse7,'Se',posse8,'Se',posse9);


tb.set_unit_cell('Bi',pos_bi1,'Bi',pos_bi2,'Se',pos_se1,'Se',pos_se3,'Se',pos_se3); %give unit cell atoms and their locations



tb.set_orbital('px,py,pz');
tb.spin_orbit_coupling('true');
tb.set_metric_unit('A');
tb.add_hrfile('bi2se3_hr.dat');
tb.set_fermi_level(4.4195);



% G 0.00000 0.00000 0.0000 Z 0.00000 0.00000 0.5000
% Z 0.00000 0.00000 0.5000 F 0.50000 0.50000 0.0000
% F 0.50000 0.50000 0.0000 G 0.00000 0.00000 0.0000
% G 0.00000 0.00000 0.0000 L 0.50000 0.00000 0.0000

if(0)
G = {[0 0 0],'$$\Gamma$$'};
Z = {2*pi*[0 0 0.5],'Z'};
F = {2*pi*[0.5 0.5 0],'F'};
L = {2*pi*[0.5 0 0],'L'};

K = {[0.33 0.67 0],'K'};
M = {[0.5 0.5 0],'M'};


precision = 150;

fig_hsym = figure("Name","High Symmetry Points Figure");
plts = tb.plot_high_symmetry_points(fig_hsym,precision,G,Z,F,G,L);

%plts = tb.plot_high_symmetry_points(fig_hsym,precision,K,G,M);
end

if(0)

range = pi; 
precision = 3;
k = tb.set_kvector(-range,range,precision);
%fig_band = figure("Name","Energy Band Figure");
%surfaces = tb.plot_energy_band(fig_band,k,'surface','EdgeColor','None');

fig_fs = figure("Name","Fermi Surface");
fs = tb.plot_fermi_surface(fig_fs,k,30);

end


if(0)
range = pi; 
precision = 20;
k = tb.set_kvector(-range,range,precision);
fig_dos = figure("Name","Density of States");
ce = tb.plot_dos(fig_dos,k);
end



if(1)
lat_f = figure("Name","Lattice Figure");
gp = lattice_drawer(lat_f,20,20,20);
gp.xaxis_symmetric();
gp.yaxis_symmetric();
gp.zaxis_symmetric();


%Gallium covalent radius = 1.26 ang, 1.36 atomic
%As covalent radius = 1.19 ang, 1.14 atomic

atoma = gp.draw('sphere rgb:ff8000',0,0,0,1.36/2,'Visible','off');
atomb = gp.draw('sphere rgb:ff007f',0,0,0,1.14/2,'Visible','off');
bond = gp.draw('line black',0,0,0,0,'Visible','off','LineWidth',8,'LineStyle','-');
%rgb:FF55FF
bonds = {bond};
atoms = {atoma,atomb};
from_to = -5:7;
tb.plot_lattice(gp,"bonds",bonds,"atoms",atoms,"x",from_to,"y",from_to,"z",from_to);


legend({'Bi','Se'},'Location','southwest','Interpreter','Latex','FontSize',24) 

%an = a;
%gp.draw('cuboid',an,an,an,an,an,an,'FaceColor','None','LineWidth',1);

end


