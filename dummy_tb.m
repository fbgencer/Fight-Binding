clear all;
close all;
clc;
a = 1e-10;

%primitive vectors
a1 = (a/2).*[0, 1, 1];
a2 = (a/2).*[1, 0, 1];
a3 = (a/2).*[1, 1, 0];


tb = tightbinding(3,a1,a2,a3);% Start with name and primitive vectors

tb.set_unit_cell('a',[0 0 0],'c',[1 1 1].*(a/4) );

tb.set_orbital('3s','p');
tb.spin_orbit_coupling('true');

%Diagonal terms
tb.add_hopping(10 ,'a','c',[0 0 0],'3s','3s','sym',sym('E_sa','real'));
tb.add_soc(1i*3,'a-','c-','3s','3s','sym','G');


tb.symbolic_hamiltonian('exact')