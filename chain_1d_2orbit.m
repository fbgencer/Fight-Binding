% 1 dimensional chain with two orbitals
clear; clc; close all;
%Chain consts
Eo = 1;
t = 2;
a0 = 1.42*1e-10;
a = a0;

%primitive vectors
a1 = [a, 0, 0];

%We have two orbitals in one unit cell
tb = tightbinding('1D-Chain',2,[a1]);

tb.add_hopping(Eo,1,1,[0]);
tb.add_hopping(Eo,2,2,[0]);
tb.add_hopping(-t,1,2,[0]);
tb.add_hopping(-t,2,1,[0]);

tb.add_hopping(-t,2,1,[1]);
tb.add_hopping(-t,1,2,[-1]);


range = 2*pi/a;
len = 80;
tb.set_kvector(-range,range,len);
tb.calculate_band();
f = tb.plot_band();