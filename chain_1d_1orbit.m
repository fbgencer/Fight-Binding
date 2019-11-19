% 1 dimensional chain with two orbitals
clear; clc; close all;
%Chain consts
Eo = 0;
t = 0.5;
a0 = 1.42*1e-10;
a = a0;

%primitive vectors
a1 = [a, 0, 0];
%We have two orbitals in one unit cell
% o-o-o-o- : chain

tb = tightbinding(1,a1);% Start with name and primitive vectors
tb.set_unit_cell('A',[0 0]); %give unit cell atoms and their locations
tb.set_orbital('s');

tb.add_hopping(Eo,'A','A',[0],'sym','Eo'); %hem A hem de 1 çalışıyor, bazı orbitalleri böyle verebiliyoruz
tb.add_hopping(-t,1,1,[1],'sym','-t');
%tb.add_hopping(-t,1,1,[-1]);

tb.symbolic_hamiltonian()

if(1)
lat_f = figure("Name","Lattice Figure");
gp = lattice_drawer(lat_f,20,20);
atoma = gp.draw('circle red',0,0,0.25,'Visible','off');
bond = gp.draw('line black',0,0,0,0,'Visible','off','LineWidth',1);
bonds = {bond};
atoms = {atoma};
tb.plot_lattice(gp,"bonds",bonds,"atoms",atoms);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(0)
range = 2*pi/a; 
precision = 100;
k = tb.set_kvector(-range,range,precision);
fig_band = figure("Name","Energy Band Figure");
surfaces = tb.plot_energy_band(fig_band,k,'surface','EdgeColor','None');

% rp = lattice_drawer(fig_band);
% lin = rp.draw('line rgb:FF8000',0,0,0,0,'Visible','off','ZData',0.5,'LineWidth',2);
% tb.plot_brillouin_zone(rp,'plot points','plot lines',lin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fig_hsym = figure("Name","High Symmetry Points Figure");
f2 = tb.plot_high_symmetry_points(fig_hsym,[0 0 0],[3.2/a 0 0]);
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