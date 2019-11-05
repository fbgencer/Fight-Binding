clear; clc; close all;
%graphene consts
Eo = 1;
t = 2;
a0 = 1.42*1e-10;
a = a0;

%primitive vectors
a1 = [a, 0, 0];
a2 = [0, a, 0];

tb = tightbinding('Square Lattice',a1,a2);
tb.set_unit_cell('A',[0 0],'B',[a/2 a/2]); %give unit cell atoms and their locations


tb.add_hopping(Eo,1,1,[0 0]);  % 0 to 0
tb.add_hopping(Eo,2,2,[0 0]);  % 0 to 0

tb.add_hopping(-t,1,2,[0 0]);
tb.add_hopping(-t,1,2,[-1 0]);
tb.add_hopping(-t,1,2,[0 -1]);
tb.add_hopping(-t,1,2,[-1 -1]);

tb.add_hopping(-t,2,1,[0 0]);
tb.add_hopping(-t,2,1,[1 0]);
tb.add_hopping(-t,2,1,[0 1]);
tb.add_hopping(-t,2,1,[1 1]);



lat_f = figure(3);
gp = lattice_drawer(lat_f,10,10);
atoma = gp.draw('circle blue',0,0,0.15,'Visible','off');
atomb = gp.draw('circle magenta',0,0,0.25,'Visible','off');
bond = gp.draw('line red',0,0,0,0,'Visible','off');
type_struct.bonds = {bond};
type_struct.atoms = {atoma,atomb};
tb.plot_lattice(gp,type_struct);


range = pi/a;
len = 80;
tb.set_kvector(-range,range,len);
tb.calculate_band();
fig_band = figure(1);
tb.plot_band(fig_band);
