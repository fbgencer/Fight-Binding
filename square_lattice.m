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
tb.set_unit_cell('A',[0 0]); %give unit cell atoms and their locations


tb.add_hopping(Eo,1,1,[0 0]);  % 0 to 0

tb.add_hopping(-t,1,1,[1 0]);  % 0 to 1 with -a1 vector
tb.add_hopping(-t,1,1,[-1 0]);   % 0 to 2 with a1 vector
tb.add_hopping(-t,1,1,[0 1]);  % 0 to 3
tb.add_hopping(-t,1,1,[0 -1]);  % 0 to 4


lat_f = figure(3);
gp = lattice_drawer(lat_f,10,10);
atoma = gp.draw('circle blue',0,0,0.25,'Visible','off');
bond = gp.draw('line red',0,0,0,0,'Visible','off');
type_struct.bonds = {bond};
type_struct.atoms = {atoma};
tb.plot_lattice(gp,0,0,type_struct);


range = 2*pi/a;
len = 80;
tb.set_kvector(-range,range,len);
tb.calculate_band();
plot_f = figure(1);
tb.plot_band(plot_f);

