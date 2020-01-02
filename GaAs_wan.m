% GaAs



clear; clc; close all;

a = 1;

%primitive vectors
a1 = (a/2).*[0, 1, 1];
a2 = (a/2).*[1, 0, 1];
a3 = (a/2).*[1, 1, 0];


tb = tightbinding(3,a1,a2,a3);% Start with name and primitive vectors

tb.set_unit_cell('a',[0 0 0],'c',[1 1 1].*(a/4) );
tb.set_orbital('3s','px,py,pz','1,2,3,4,5');
tb.set_metric_unit('m');
tb.spin_orbit_coupling('false');
tb.add_hrfile('GaAs_hr.dat');
tb.set_fermi_level(6);

if(0)
lat_f = figure("Name","Lattice Figure");
gp = lattice_drawer(lat_f,3,3,3);


atoma = gp.draw('sphere rgb:33FFFF',0,0,0,0.1,'Visible','off');
atomb = gp.draw('sphere black',0,0,0,0.1,'Visible','off');
bond = gp.draw('line yellow',0,0,0,0,'Visible','off','LineWidth',4,'LineStyle','-');
%rgb:FF55FF
bonds = {bond};
atoms = {atoma,atomb};
from_to = 1:2;
tb.plot_lattice(gp,"bonds",bonds,"atoms",atoms,"x",from_to,"y",from_to,"z",from_to);


an = a;
gp.draw('cuboid',an,an,an,an,an,an,'FaceColor','None','LineWidth',1);

end


if(0)
% begin kpoint_path
% L   0.50000     0.50000     0.50000   G   0.00000     0.00000     0.00000
% G   0.00000     0.00000     0.00000   X   0.50000     0.00000     0.50000
% X   0.50000     0.00000     0.50000   W   0.50000     0.25000     0.75000
% W   0.50000     0.25000     0.75000   L   0.50000     0.50000     0.50000
% L   0.50000     0.50000     0.50000   K   0.75000     0.37500     0.37500
% K   0.75000     0.37500     0.37500   G   0.00000     0.00000     0.00000
% end kpoint_path


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


if(1)
range = pi; 
precision = 50;
k = tb.set_kvector(-range,range,precision);
fig_dos = figure("Name","Density of States");
ce = tb.plot_dos(fig_dos,k);
end
