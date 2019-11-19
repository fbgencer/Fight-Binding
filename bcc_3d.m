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


tb = tightbinding(3,a1,a2,a3);% Start with name and primitive vectors

tb.set_unit_cell('A',[0 0 0]);

tb.add_hopping(Eo,1,1,[0 0 0]);

tb.add_hopping(-t,1,1,[1 0 0]);
tb.add_hopping(-t,1,1,[0 1 0]);
tb.add_hopping(-t,1,1,[0 0 1]);
tb.add_hopping(-t,1,1,[1 1 1]);

%tb.add_hopping(-t,1,1,[-1 0 0]);
%tb.add_hopping(-t,1,1,[0 -1 0]);
%tb.add_hopping(-t,1,1,[0 0 -1]);
%tb.add_hopping(-t,1,1,[-1 -1 -1]);

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


if(1)
lat_f = figure("Name","Lattice Figure");
gp = lattice_drawer(lat_f,6,6,6);
atoma = gp.draw('point red',0,0,0.25,'Visible','off');
bond = gp.draw('line black',0,0,0,0,'Visible','off');
bonds = {bond};
atoms = {atoma};
tb.plot_lattice(gp,"bonds",bonds,"atoms",atoms,"x",-3:3,"y",-1:1,"z",-1:1);
end


if(0)

range = 2*pi/a; 
precision = 50;
k = tb.set_kvector(-range,range,precision);
fig_band = figure("Name","Energy Band Figure");
surfaces = tb.plot_energy_band(fig_band,k,'surface','EdgeColor','None');

colorbar;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(1)
Gamma = [0 0 0];
P = [pi/a pi/a pi/a ];
N = [0 pi/a pi/a];
H = [0 0 2*pi/a];


fig_hsym = figure("Name","High Symmetry Points Figure");
tb.plot_high_symmetry_points(fig_hsym,P,Gamma,N,H,Gamma);
xlabel(fig_hsym.CurrentAxes,'$$P \Gamma N H \Gamma $$','Interpreter','Latex');
xticks(fig_hsym.CurrentAxes,[0 80 160 240 320]);
xticklabels(fig_hsym.CurrentAxes,{'P','\Gamma','N','H','\Gamma'});

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