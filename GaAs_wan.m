% GaAs



clear; clc; close all;

a = 5.65; % unit of arm
a1 = (a/2).*[0, 1, 1];
a2 = (a/2).*[1, 0, 1];
a3 = (a/2).*[1, 1, 0];


tb = tightbinding(3,a1,a2,a3);% Start with name and primitive vectors

tb.set_atom_types('Ga','As');
tb.set_unit_cell('Ga',[0 0 0],'As',[0.25 0.25 0.25].*a );
tb.set_orbital('1s','2s','2px,2py,2pz','3s','3px,3py,3pz');
tb.set_metric_unit('A');
tb.spin_orbit_coupling('true');
tb.add_hrfile('GaAs_hr.dat');
tb.set_fermi_level(7.1912);

if(0)
lat_f = figure("Name","Lattice Figure");
gp = lattice_drawer(lat_f,15,15,15);

%Gallium covalent radius = 1.26 ang, 1.36 atomic
%As covalent radius = 1.19 ang, 1.14 atomic

atoma = gp.draw('sphere rgb:4522e4',0,0,0,1.36/2,'Visible','off');
atomb = gp.draw('sphere rgb:66ff66',0,0,0,1.14/2,'Visible','off');
bond = gp.draw('line rgb:ff8000',0,0,0,0,'Visible','off','LineWidth',8,'LineStyle','-');
%rgb:FF55FF
bonds = {bond};
atoms = {atoma,atomb};
from_to = 5:12;
tb.plot_lattice(gp,"bonds",bonds,"atoms",atoms,"x",from_to,"y",from_to,"z",from_to,'unit vector',[11.3 11.3 11.3]);


an = a;
gp.draw('cuboid',an,an,an,an,an,an,'FaceColor','None','LineWidth',1);

legend({'Ga','As'},'Location','southwest','Interpreter','Latex','FontSize',24)
end 

if(1)

L = {2*pi*[0.5 0.5 0.5],'L'};
G = {[0 0 0],'$$\Gamma$$'};
X = {2*pi*[0.5 0 0.5],'X'};
W = {2*pi*[0.5 0.25 0.75],'W'};
K = {2*pi*[0.75 0.375 0.375],'K'};

precision = 30;

fig_hsym = figure("Name","High Symmetry Points Figure");
tb.plot_high_symmetry_points(fig_hsym,precision,L,G,X);
%,W,L,K,G
end


if(0)
range = pi; 
precision = 10;
k = tb.set_kvector(-range,range,precision);
fig_dos = figure("Name","Density of States");
ce = tb.plot_dos(fig_dos,k,'k -','LineWidth',1.5);
end



if(0)

range = pi; 
precision = 5;
k = tb.set_kvector(-range,range,precision);
%fig_band = figure("Name","Energy Band Figure");
%surfaces = tb.plot_energy_band(fig_band,k,'surface','EdgeColor','None');

%fig_fs = figure("Name","Fermi Surface");
for i = 22:36
f = figure(i);
fs = tb.plot_fermi_surface(f,k,i);
end
end