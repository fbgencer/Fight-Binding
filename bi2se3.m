clear all;
clc; close all;

a1 = [-2.069  -3.583614  0.000000];
a2 = [2.069  -3.583614  0.000000];
a3 = [0.000   2.389075  9.546667];

tb = tightbinding(3,a1,a2,a3);% Start with dimension and primitive vectors

pos_bi1 = [0.000000    -1.194539    6.654027];
pos_bi2 = [0.000000    -3.583614    2.892640];
pos_se1 = [0.000000     1.194538    4.773334];
pos_se2 = [0.000000    -1.194538    1.126507];
pos_se3 = [0.000000    -3.583615    8.420160];


tb.set_atom_types('Bi','Se');
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


if(1)
lat_f = figure("Name","Lattice Figure");
gp = lattice_drawer(lat_f,20,20,20);
gp.axis_symmetric();

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




if(1)
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


if(1)
range = pi; 
precision = 10;
k = tb.set_kvector(-range,range,precision);
fig_dos = figure("Name","Density of States");
ce = tb.plot_dos(fig_dos,k,'k -','LineWidth',1.5);
end




