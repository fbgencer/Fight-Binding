% NbSe2
clear all; clc; close all;

%primitive vectors
a1 = [2.9439482089090570 -1.6996892908939607 0.0000000000000000];
a2 = [0.0000000000000000  3.3993785817879214   0.0000000000000000];
a3 = [0.0000000000000000  0.0000000000000000 24.0320930230956016];


tb = tightbinding(3,a1,a2,a3);% Start with name and primitive vectors
Nb1=[0.000000    0.000000    12.016047];
Se1=[0.981316    1.699689    13.682898];
Se2=[0.981316    1.699689    10.349196];

tb.set_atom_types('Nb','Se');

tb.set_unit_cell('Nb',Nb1,'Se',Se1,'Se',Se2);

%tb.set_orbital('Nb:d1,d2,d3,d4,d5','Se:px,py,pz');
tb.set_matrix_size(11);
tb.set_metric_unit('A');
tb.spin_orbit_coupling('false');
tb.add_hrfile('NbSe2_hr_3site.dat');
tb.set_fermi_level(-1.1865);



 
if(1)

G = {[0 0 0],'$$\Gamma$$'};
M = {2*pi*[0.5 0 0],'M'};
K = {2*pi*[0.33 0.33 0],'K'};
L = {2*pi*[0.5 0 0.5],'L'};

precision = [100,30,100];

fig_hsym = figure("Name","High Symmetry Points Figure");
tb.plot_high_symmetry_points(fig_hsym,precision,G,M,K,G);
ax = fig_hsym.CurrentAxes;
ax.YLim = [-6 5];
end


if(1)
range = pi; 
precision = 10;
k = tb.set_kvector(-range,range,precision);
fig_dos = figure("Name","Density of States");
ce = tb.plot_dos(fig_dos,k,'k -','LineWidth',1.5);
end

if(1)
lat_f = figure("Name","Lattice Figure");
gp = lattice_drawer(lat_f,20,20,30);
gp.axis_symmetric();

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
%,"z",-6:6,"x",-3:3,"y",-3:3
legend({'Nb','Se'},'Location','southwest','Interpreter','Latex','FontSize',24) 

end