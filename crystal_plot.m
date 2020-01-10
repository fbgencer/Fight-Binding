clear all; close all; clc;
a = 5.65; %arm
a1 = (a/2).*[0, 1, 1];
a2 = (a/2).*[1, 0, 1];
a3 = (a/2).*[1, 1, 0];

tb = tightbinding(3,a1,a2,a3);% Start with name and primitive vectors

tb.set_unit_cell('a',[0 0 0],'c',[0.25 0.25 0.25].*a );
tb.set_orbital('3s','px,py,pz','4s');
tb.set_metric_unit('A');

lat_f = figure("Name","Lattice Figure");
gp = lattice_drawer(lat_f,15,15,15);


%Gallium covalent radius = 1.26 ang, 1.36 atomic
%As covalent radius = 1.19 ang, 1.14 atomic

atoma = gp.draw('sphere rgb:4589e4',0,0,0,1.36/2,'Visible','off');
atomb = gp.draw('sphere black',0,0,0,1.14/2,'Visible','off');
bond = gp.draw('line red',0,0,0,0,'Visible','off','LineWidth',4,'LineStyle','-');
%rgb:FF55FF
bonds = {bond};
atoms = {atoma,atomb};
from_to = 5:12;
tb.plot_lattice(gp,"bonds",bonds,"atoms",atoms,"x",from_to,"y",from_to,"z",from_to);


an = a;
gp.draw('cuboid',an,an,an,an,an,an,'FaceColor','None','LineWidth',1);
