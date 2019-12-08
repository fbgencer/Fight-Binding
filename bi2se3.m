clear all;
clc; close all;

a1 = [-2.069  -3.583614  0.000000];
a2 = [2.069  -3.583614  0.000000];
a3 = [0.000   2.389075  9.546667];

tb = tightbinding(3,a1,a2,a3);% Start with dimension and primitive vectors

pos_bi1 = [0.3990    0.3990    0.6970];
pos_bi2 = [0.6010    0.6010    0.3030];
pos_se1 = [0.0000    0.0000    0.5000];
pos_se2 = [0.2060    0.2060    0.1180];
pos_se3 = [0.7940    0.7940    0.8820];


tb.set_unit_cell('Bi1',pos_bi1,'Bi2',pos_bi2,'Se1',pos_se1,'Se2',pos_se3,'Se3',pos_se3); %give unit cell atoms and their locations
tb.set_orbital('px,py,pz');
tb.spin_orbit_coupling('true');
tb.add_hrfile('bi2se3_hr.dat');
tb.set_fermi_level(4.4195);


% G 0.00000 0.00000 0.0000 Z 0.00000 0.00000 0.5000
% Z 0.00000 0.00000 0.5000 F 0.50000 0.50000 0.0000
% F 0.50000 0.50000 0.0000 G 0.00000 0.00000 0.0000
% G 0.00000 0.00000 0.0000 L 0.50000 0.00000 0.0000

if(0)
G = [0 0 0];
Z = [0 0 0.5];
F = [0.5 0.5 0];
L = [0.5 0 0];

precision = 71;

fig_hsym = figure("Name","High Symmetry Points Figure");
plts = tb.plot_high_symmetry_points(fig_hsym,precision,G,Z,F,G,L);

xlabel(fig_hsym.CurrentAxes,'$$G Z F G L$$','Interpreter','Latex');

xticks(fig_hsym.CurrentAxes,0:precision:5*precision);
xticklabels(fig_hsym.CurrentAxes,{'G','Z','F','G','L'});
ylabel(fig_hsym.CurrentAxes,'$$Energy(eV)$$','Interpreter','Latex')
grid();
%fig_hsym.CurrentAxes.YLim = [-3 3];
%for ix = 1:numel(plts), plts{ix}.Color = 'red'; end
end


if(1)
range = 2*pi; 
precision = 21;
k = tb.set_kvector(-range,range,precision);
fig_dos = figure("Name","Density of States");
ce = tb.plot_dos(fig_dos,k);
end



if(0)
lat_f = figure("Name","Lattice Figure");
gp = lattice_drawer(lat_f,3,3,3);
gp.xaxis_symmetric();
gp.yaxis_symmetric();
gp.zaxis_symmetric();

atom1 = gp.draw('point blue',0,0,0.1,'Visible','off');
atom2 = gp.draw('point black',0,0,0.1,'Visible','off');
bond = gp.draw('line rgb:FF55FF',0,0,0,0,'Visible','off','LineWidth',0.1,'LineStyle','--');
bonds = {bond};
atoms = {atom1,atom1,atom2,atom2,atom2};
from_to = -2:2;
tb.plot_lattice(gp,"bonds",bonds,"atoms",atoms,"x",from_to,"y",from_to,"z",from_to);
gp.set_xlabel('X');
gp.set_ylabel('Y');
gp.set_zlabel('Z');
grid;

an = a;
gp.draw('cuboid',an,an,an,an,an,an,'FaceColor','None','LineWidth',1);
end


