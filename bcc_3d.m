% 1 dimensional chain with two orbitals
clear; clc; close all;
%Chain consts
Eo = 0;
t = 1;
a0 = 1.42*1e-10;
a = a0;

%primitive vectors
a1 = (a/2).*[1, 1, -1];
a2 = (a/2).*[-1, 1, 1];
a3 = (a/2).*[1, -1, 1];


tb = tightbinding('Bcc',a1,a2,a3);% Start with name and primitive vectors


tb.set_unit_cell('A',[0 0 0]);
tb.add_hopping(Eo,1,1,[0 0 0]);
tb.add_hopping(-t,1,1,[1 0 0]);
tb.add_hopping(-t,1,1,[0 1 0]);
tb.add_hopping(-t,1,1,[0 0 1]);
tb.add_hopping(-t,1,1,[-1 0 0]);
tb.add_hopping(-t,1,1,[0 -1 0]);
tb.add_hopping(-t,1,1,[0 0 -1]);
tb.add_hopping(-t,1,1,[1 1 1]);
tb.add_hopping(-t,1,1,[-1 -1 -1]);

%Bu düzgün çalışmadı
%tb.set_unit_cell('1',[0,0,0],'2',[-a/2, a/2, -a/2],'3',[a/2, a/2, -a/2],'4',[-a/2, a/2, a/2], ...
%	'5',[a/2, a/2, a/2],'6',[-a/2, -a/2, -a/2],'7',[a/2, -a/2, -a/2],'8',[-a/2, -a/2, a/2],'9',[a/2, -a/2, a/2]); %give unit cell atoms and their locations
%tb.add_hopping(Eo,1,1,[0 0 0]);
%tb.add_hopping(-t,1,2,[0 0 0]);
%tb.add_hopping(-t,1,3,[0 0 0]);
%tb.add_hopping(-t,1,4,[0 0 0]);
%tb.add_hopping(-t,1,5,[0 0 0]);
%tb.add_hopping(-t,1,6,[0 0 0]);
%tb.add_hopping(-t,1,7,[0 0 0]);
%tb.add_hopping(-t,1,8,[0 0 0]);
%tb.add_hopping(-t,1,9,[0 0 0]);


if(0)
lat_f = figure(3);
gp = lattice_drawer(lat_f,5,5,5);
atoma = gp.draw('point red',0,0,0.25,'Visible','off');
atomb = gp.draw('point black',0,0,0.15,'Visible','off');
bond = gp.draw('line black',0,0,0,0,'Visible','off');
type_struct.bonds = {bond};
type_struct.atoms = {atoma};
%tb.plot_lattice(gp,type_struct);
tb.plot_lattice(gp,type_struct);
%tb.plot_only_atoms(gp,type_struct);
%tb.plot_only_bonds(gp,type_struct);
end


if(1)

range = 2*pi/a;
len = 50;
tb.set_kvector(-range,range,len);
tb.calculate_band();
fig_band = figure(1);
f = tb.plot_band(fig_band);
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gamma = [0 0 0];
P = [pi/a pi/a pi/a ];
N = [0 pi/a pi/a];
H = [0 0 2*pi/a];


%k_fig = figure(2);
%f2 = tb.plot_high_symmetry_points(k_fig,P,Gamma,N,H,Gamma);
%xlabel(f2.CurrentAxes,'$$P \Gamma M K$$','Interpreter','Latex');
%xticks(f2.CurrentAxes,[0 80 160 240 320]);
%xticklabels(f2.CurrentAxes,{'P','\Gamma','N','H','\Gamma'});

end

if(0)
%Hand solution
Gamma = [0 0 0];
P = [pi/a pi/a pi/a ];
N = [0 pi/a pi/a];
H = [0 0 2*pi/a];

points = {P,Gamma,N,H,Gamma};

precision = 50;
kx = [];
ky = [];
kz = [];
for i = 1:size(points,2)-1
    kx = [kx,linspace(points{i}(1),points{i+1}(1),precision)];
    ky = [ky,linspace(points{i}(2),points{i+1}(2),precision)];
    kz = [kz,linspace(points{i}(3),points{i+1}(3),precision)];
end


E = Eo - 8.*t.*cos(kx*a/2).*cos(ky*a/2).*cos(kz*a/2);

plot(E);

end