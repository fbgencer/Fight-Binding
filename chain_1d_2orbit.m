% 1 dimensional chain with two orbitals
clear; clc; close all;
%Chain consts
Eo = 0;
t = 0.5;
a0 = 1.42*1e-10;
a = a0;

%primitive vectors
a1 = [a, 0, 0];
a2 = [0 0 0];
%We have two orbitals in one unit cell
% o-O--o-O--o-O--o-O-- : chain

tb = tightbinding('1D-Chain',a1,a2);% Start with name and primitive vectors
tb.set_unit_cell('A',[-a0/4 0],'B',[a0/4 0]); %give unit cell atoms and their locations


tb.add_hopping(Eo,1,1,[0 0]);
tb.add_hopping(Eo,2,2,[0 0]);
tb.add_hopping(-t,1,2,[0 0]);
tb.add_hopping(-t,2,1,[0 0]);

tb.add_hopping(-t,2,1,[1 0]);
tb.add_hopping(-t,1,2,[-1 0]);

if(1)
lat_f = figure(3);
gp = lattice_drawer(lat_f,20,20);
atoma = gp.draw('circle blue',0,0,0.25,'Visible','off');
atomb = gp.draw('circle magenta',0,0,0.1,'Visible','off');
bond = gp.draw('line red',0,0,0,0,'Visible','off');
type_struct.bonds = {bond};
type_struct.atoms = {atoma,atomb};
%tb.plot_lattice(gp,type_struct);

a1 = a1.*1e10;
a2 = a2.*1e10;


for e1 = -1:1
	for e2 = -1:1
        cx = e1 * a1(1) + e2 * a2(1);
        cy = e1 * a1(2) + e2 * a2(2);
        
		tb.plot_only_bonds(gp,cx,cy,type_struct);
		tb.plot_only_atoms(gp,cx,cy,type_struct);
	end
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(0)
range = 2*pi/a;
len = 80;
tb.set_kvector(-range,range,len);
tb.calculate_band();
fig_band = figure(1);
tb.plot_band(fig_band);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


k_fig = figure(2);
f2 = tb.plot_high_symmetry_points(k_fig,[0 0 0],[3.2/a 0 0]);
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