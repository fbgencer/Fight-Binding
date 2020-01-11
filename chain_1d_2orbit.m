% 1 dimensional chain with two orbitals
clear; clc; close all;
%Chain consts
Eo = 0;
t = 0.5;
a0 = 1;
a = a0;

%primitive vectors
a1 = [a, 0, 0];
%We have two atoms in one unit cell
% o-O--o-O--o-O--o-O-- : chain

tb = tightbinding(1,a1);% Start with dimension and primitive vectors
tb.set_unit_cell('A',[-a0/4 0],'B',[a0/4 0]); %give unit cell atoms and their locations
%tb.set_orbital('s','p');
tb.set_metric_unit('A');

tb.add_hopping(Eo,1,1,[0]);
tb.add_hopping(Eo,2,2,[0]);
tb.add_hopping(-t,1,2,[0]);
%tb.add_hopping(-t,2,1,[0]);

tb.add_hopping(-t,2,1,[1]);
%tb.add_hopping(-t,1,2,[-1]);

if(1)
lat_f = figure("Name","Lattice Figure");
gp = lattice_drawer(lat_f,16,10);
gp.xaxis_symmetric();
gp.yaxis_symmetric();
atoma = gp.draw('circle blue',0,0,0.15,'Visible','off');
atomb = gp.draw('circle red',0,0,0.15,'Visible','off');
bond = gp.draw('line black',0,0,0,0,'Visible','off');
bonds = {bond};
atoms = {atoma,atomb};
tb.plot_lattice(gp,"bonds",bonds,"atoms",atoms,"x",-6:6,"y",-2:2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(0)
range = pi; 
precision = 300;
k = tb.set_kvector(-range,range,precision);
fig_dos = figure("Name","Density of States");
ce = tb.plot_dos(fig_dos,k);
end



if(1)
range = 2*pi/a; 
precision = 100;
k = tb.set_kvector(-range,range,precision);
%fig_band = figure("Name","Energy Band Figure");
%surfaces = tb.plot_energy_band(fig_band,k,'surface','EdgeColor','None');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig_hsym = figure("Name","High Symmetry Points Figure");
any_point = {[3.2/a 0 0],'L'};
G = {[0 0 0],'$$\Gamma$$'};
f2 = tb.plot_high_symmetry_points(fig_hsym,precision,G,any_point);

end

%Hand solution
%Use only for verification
%e1 = Eo;
%e2 = Eo;
%k = linspace(0,3.2/a,80);
%E = 0.5*((e1+e2)+sqrt( (e1-e2).^2 + 8*t*t*(1+cos(k*a)) ));
%En  = 0.5*((e1+e2)-sqrt( (e1-e2).^2 + 8*t*t*(1+cos(k*a)) ));

%figure(3);
%plot(E)
%hold on;
%plot(En);