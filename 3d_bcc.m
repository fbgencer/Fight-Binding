% 1 dimensional chain with two orbitals
clear; clc; close all;
%Chain consts
Eo = 0;
t = 0.5;
a0 = 1.42*1e-10;
a = a0;

%primitive vectors
a1 = (a/2)*[1, 1, -1];
a2 = (a/2)*[-1, 1, 1];
a3 = (a/2)*[1, -1, 1];


tb = tightbinding('Bcc',a1,a2,a3);% Start with name and primitive vectors
tb.set_unit_cell('A',[-a0/2 0],'B',[a0/2 0]); %give unit cell atoms and their locations


tb.add_hopping(Eo,1,1,[0 0]);
tb.add_hopping(Eo,2,2,[0 0]);
tb.add_hopping(-t,1,2,[0 0]);
tb.add_hopping(-t,2,1,[0 0]);

tb.add_hopping(-t,2,1,[1 0 ]);
tb.add_hopping(-t,1,2,[-1 0]);



range = 2*pi/a;
len = 80;
tb.set_kvector(-range,range,len);
tb.calculate_band();
fig_band = figure(1);
tb.plot_band(fig_band);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k_fig = figure(2);
f2 = tb.plot_high_symmetry_points(k_fig,[0 0 0],[3.2/a 0 0]);

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