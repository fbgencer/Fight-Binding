clear; clc; close all;
%graphene consts
Eo = 1;
t = 2;
a0 = 1.42*1e-10;
a = a0;

%primitive vectors
a1 = [a, 0, 0];
a2 = [0, a, 0];

tb = tightbinding(2,a1,a2);
tb.set_unit_cell('A',[0 0]); %give unit cell atoms and their locations


tb.add_hopping(Eo,1,1,[0 0]);  % 0 to 0

tb.add_hopping(-t,1,1,[1 0]);  % 0 to 1 with -a1 vector
tb.add_hopping(-t,1,1,[-1 0]);   % 0 to 2 with a1 vector
tb.add_hopping(-t,1,1,[0 1]);  % 0 to 3
tb.add_hopping(-t,1,1,[0 -1]);  % 0 to 4


fig_lat = figure("Name","Lattice Figure");
gp = lattice_drawer(fig_lat,10,10);
atoma = gp.draw('crect blue',0,0,0.3,0.3,'Visible','off');
bond = gp.draw('line red',0,0,0,0,'Visible','off');
bonds = {bond};
atoms = {atoma};
tb.plot_lattice(gp,"bonds",bonds,"atoms",atoms);


range = 2*pi/a;
len = 100;
k = tb.set_kvector(-range,range,len);
fig_band= figure("Name","Energy Band Figure");
surfaces = tb.plot_energy_band(fig_band,k,'surface','EdgeColor','None');

rp = lattice_drawer(fig_band);
pt = rp.draw('point rgb:660066',0,0,10,'Visible','off');
lin = rp.draw('line rgb:FF8000',0,0,0,0,'Visible','off','ZData',10,'LineWidth',2);
tb.plot_brillouin_zone(rp,'plot points',pt,'plot lines',lin,'plot coordinates');