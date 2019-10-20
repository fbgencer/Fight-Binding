clear; clc; close all;
%graphene consts
Eo = 1;
t = 2;
a0 = 1.42*1e-10;
a = a0;

%primitive vectors
a1 = [a, 0, 0];
a2 = [a/2, a/2, 0];

tb = tightbinding('Square Lattice 2 Atom',2,[a1;a2]);

tb.add_hopping(Eo,1,1,[0 0]);  % 0 to 0
tb.add_hopping(Eo,2,2,[0 0]);  % 0 to 0

tb.add_hopping(-t,1,1,[1 0]);  % 0 to 1 with -a1 vector
tb.add_hopping(-t,1,1,[-1 0]);   % 0 to 2 with a1 vector
tb.add_hopping(-t,1,1,[0 1]);  % 0 to 3
tb.add_hopping(-t,1,1,[0 -1]);  % 0 to 4



range = 2*pi/a;
len = 80;
tb.set_kvector(-range,range,len);
tb.calculate_band();
f = tb.plot_band();

