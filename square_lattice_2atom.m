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

fig_lat = figure("Name","Lattice Figure");
gp = lattice_drawer(fig_lat,10,10);
atoma = gp.draw('circle rgb:660066',0,0,0.15,'Visible','off');
atomb = gp.draw('circle rgb:00944c',0,0,0.25,'Visible','off');
bond = gp.draw('line k',0,0,0,0,'Visible','off');
bonds = {bond};
atoms = {atoma,atomb};
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