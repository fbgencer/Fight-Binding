% Gaas

clear; clc; close all;

%Chain consts
Esa=-8.3431;
Epa=1.0414;
Esc=-2.6569;
Epc=3.6686;
Esea=8.5914;
Esec=6.7386;
Ess=-6.4513;
Exx=1.9546;
Exy=5.0779;
Esapc=4.4800;
Epasc=5.7839;
Eseapc=4.8422;

a = 1.42*1e-10; 

%primitive vectors
a1 = (a/2).*[0, 1, 1];
a2 = (a/2).*[1, 0, 1];
a3 = (a/2).*[1, 1, 0];


tb = tightbinding(3,a1,a2,a3);% Start with name and primitive vectors

tb.set_unit_cell('a',[0 0 0],'c',[1 1 1].*a/4);
tb.set_orbital('s','px','py','pz','s*');

%Diagonal terms
tb.add_hopping(Esa,'a','a',[0 0 0],'s','s','sym','E_sa');
tb.add_hopping(Esc,'c','c',[0 0 0],'s','s','sym','E_sc');
tb.add_hopping(Epa,'a','a',[0 0 0],'px','px','sym','E_pa');
tb.add_hopping(Epa,'a','a',[0 0 0],'py','py','sym','E_pa');
tb.add_hopping(Epa,'a','a',[0 0 0],'pz','pz','sym','E_pa');
tb.add_hopping(Epc,'c','c',[0 0 0],'px','px','sym','E_pc');
tb.add_hopping(Epc,'c','c',[0 0 0],'py','py','sym','E_pc');
tb.add_hopping(Epc,'c','c',[0 0 0],'pz','pz','sym','E_pc');
tb.add_hopping(Esea,'a','a',[0 0 0],'s*','s*','sym','E_s*a');
tb.add_hopping(Esec,'c','c',[0 0 0],'s*','s*','sym','E_s*c');

%%first row
tb.add_hopping(Ess,'a','c',[0 0 0],'s','s','sym','E_ss');
tb.add_hopping(Ess,'a','c',[1 0 0],'s','s','sym','E_ss');
tb.add_hopping(Ess,'a','c',[0 1 0],'s','s','sym','E_ss');
tb.add_hopping(Ess,'a','c',[0 0 1],'s','s','sym','E_ss');

tb.add_hopping(0,'a','a',[0 0 0],'s','px','sym','0');
tb.add_hopping(0,'a','a',[0 0 0],'s','py','sym','0');
tb.add_hopping(0,'a','a',[0 0 0],'s','pz','sym','0');

tb.add_hopping(Esapc,'a','c',[0 0 0],'s','px','sym','E_sapc');
tb.add_hopping(Esapc,'a','c',[1 0 0],'s','px','sym','E_sapc');
tb.add_hopping(-Esapc,'a','c',[0 1 0],'s','px','sym','E_sapc');
tb.add_hopping(-Esapc,'a','c',[0 0 1],'s','px','sym','E_sapc');

tb.add_hopping(Esapc,'a','c',[0 0 0],'s','py','sym','E_sapc');
tb.add_hopping(-Esapc,'a','c',[1 0 0],'s','py','sym','E_sapc');
tb.add_hopping(+Esapc,'a','c',[0 1 0],'s','py','sym','E_sapc');
tb.add_hopping(-Esapc,'a','c',[0 0 1],'s','py','sym','E_sapc');

tb.add_hopping(Esapc,'a','c',[0 0 0],'s','pz','sym','E_sapc');
tb.add_hopping(-Esapc,'a','c',[1 0 0],'s','pz','sym','E_sapc');
tb.add_hopping(-Esapc,'a','c',[0 1 0],'s','pz','sym','E_sapc');
tb.add_hopping(Esapc,'a','c',[0 0 1],'s','pz','sym','E_sapc');

tb.add_hopping(0,'a','a',[0 0 0],'s','s*','sym','0');
tb.add_hopping(0,'a','c',[0 0 0],'s','s*','sym','0');

%%2nd row
tb.add_hopping(Epasc,'c','a',[0 0 0],'s','px','sym','E_pasc');
tb.add_hopping(Epasc,'c','a',[-1 0 0],'s','px','sym','E_pasc');
tb.add_hopping(-Epasc,'c','a',[0 -1 0],'s','px','sym','E_pasc');
tb.add_hopping(-Epasc,'c','a',[0 0 -1],'s','px','sym','E_pasc');

tb.add_hopping(Epasc,'c','a',[0 0 0],'s','py','sym','E_pasc');
tb.add_hopping(-Epasc,'c','a',[-1 0 0],'s','py','sym','E_pasc');
tb.add_hopping(+Epasc,'c','a',[0 -1 0],'s','py','sym','E_pasc');
tb.add_hopping(-Epasc,'c','a',[0 0 -1],'s','py','sym','E_pasc');

tb.add_hopping(Epasc,'c','a',[0 0 0],'s','pz','sym','E_pasc');
tb.add_hopping(-Epasc,'c','a',[-1 0 0],'s','pz','sym','E_pasc');
tb.add_hopping(-Epasc,'c','a',[0 -1 0],'s','pz','sym','E_pasc');
tb.add_hopping(Epasc,'c','a',[0 0 -1],'s','pz','sym','E_pasc');

tb.add_hopping(0,'c','c',[0 0 0],'s','px','sym','0');
tb.add_hopping(0,'c','c',[0 0 0],'s','py','sym','0');
tb.add_hopping(0,'c','c',[0 0 0],'s','pz','sym','0');

tb.add_hopping(0,'c','a',[0 0 0],'s','s*','sym','0');
tb.add_hopping(0,'c','c',[0 0 0],'s','s*','sym','0');

%3rd row
tb.add_hopping(0,'a','a',[0 0 0],'px','py','sym','0');
tb.add_hopping(0,'a','a',[0 0 0],'px','pz','sym','0');

tb.add_hopping(Exx,'a','c',[0 0 0],'px','px','sym','E_xx');
tb.add_hopping(Exx,'a','c',[1 0 0],'px','px','sym','E_xx');
tb.add_hopping(Exx,'a','c',[0 1 0],'px','px','sym','E_xx');
tb.add_hopping(Exx,'a','c',[0 0 1],'px','px','sym','E_xx');

tb.add_hopping(Exy,'a','c',[0 0 0],'px','py','sym','E_xy');
tb.add_hopping(-Exy,'a','c',[1 0 0],'px','py','sym','E_xy');
tb.add_hopping(-Exy,'a','c',[0 1 0],'px','py','sym','E_xy');
tb.add_hopping(Exy,'a','c',[0 0 1],'px','py','sym','E_xy');

tb.add_hopping(Exy,'a','c',[0 0 0],'px','pz','sym','E_xy');
tb.add_hopping(-Exy,'a','c',[1 0 0],'px','pz','sym','E_xy');
tb.add_hopping(Exy,'a','c',[0 1 0],'px','pz','sym','E_xy');
tb.add_hopping(-Exy,'a','c',[0 0 1],'px','pz','sym','E_xy');

tb.add_hopping(0,'a','a',[0 0 0],'px','s*','sym','0');

tb.add_hopping(Epasc,'a','c',[0 0 0],'px','s*','sym','E_pasc');
tb.add_hopping(Epasc,'a','c',[1 0 0],'px','s*','sym','E_pasc');
tb.add_hopping(Epasc,'a','c',[0 1 0],'px','s*','sym','E_pasc');
tb.add_hopping(Epasc,'a','c',[0 0 1],'px','s*','sym','E_pasc');

%4th row
tb.add_hopping(0,'a','a',[0 0 0],'py','pz','sym','0');

tb.add_hopping(Exy,'a','c',[0 0 0],'py','px','sym','E_xy');
tb.add_hopping(-Exy,'a','c',[1 0 0],'py','px','sym','E_xy');
tb.add_hopping(-Exy,'a','c',[0 1 0],'py','px','sym','E_xy');
tb.add_hopping(Exy,'a','c',[0 0 1],'py','px','sym','E_xy');

tb.add_hopping(Exx,'a','c',[0 0 0],'py','py','sym','E_xx');
tb.add_hopping(Exx,'a','c',[1 0 0],'py','py','sym','E_xx');
tb.add_hopping(Exx,'a','c',[0 1 0],'py','py','sym','E_xx');
tb.add_hopping(Exx,'a','c',[0 0 1],'py','py','sym','E_xx');

tb.add_hopping(Exy,'a','c',[0 0 0],'py','pz','sym','E_xy');
tb.add_hopping(Exy,'a','c',[1 0 0],'py','pz','sym','E_xy');
tb.add_hopping(-Exy,'a','c',[0 1 0],'py','pz','sym','E_xy');
tb.add_hopping(-Exy,'a','c',[0 0 1],'py','pz','sym','E_xy');

tb.add_hopping(0,'a','a',[0 0 0],'py','s*','sym','0');

tb.add_hopping(Epasc,'a','c',[0 0 0],'py','s*','sym','E_pasc');
tb.add_hopping(-Epasc,'a','c',[1 0 0],'py','s*','sym','E_pasc');
tb.add_hopping(Epasc,'a','c',[0 1 0],'py','s*','sym','E_pasc');
tb.add_hopping(-Epasc,'a','c',[0 0 1],'py','s*','sym','E_pasc');

%5th row
tb.add_hopping(Exy,'a','c',[0 0 0],'pz','px','sym','E_xy');
tb.add_hopping(-Exy,'a','c',[1 0 0],'pz','px','sym','E_xy');
tb.add_hopping(Exy,'a','c',[0 1 0],'pz','px','sym','E_xy');
tb.add_hopping(-Exy,'a','c',[0 0 1],'pz','px','sym','E_xy');

tb.add_hopping(Exy,'a','c',[0 0 0],'pz','py','sym','E_xy');
tb.add_hopping(Exy,'a','c',[1 0 0],'pz','py','sym','E_xy');
tb.add_hopping(-Exy,'a','c',[0 1 0],'pz','py','sym','E_xy');
tb.add_hopping(-Exy,'a','c',[0 0 1],'pz','py','sym','E_xy');

tb.add_hopping(Exx,'a','c',[0 0 0],'pz','pz','sym','E_xx');
tb.add_hopping(Exx,'a','c',[1 0 0],'pz','pz','sym','E_xx');
tb.add_hopping(Exx,'a','c',[0 1 0],'pz','pz','sym','E_xx');
tb.add_hopping(Exx,'a','c',[0 0 1],'pz','pz','sym','E_xx');

tb.add_hopping(0,'a','a',[0 0 0],'pz','s*','sym','0');

tb.add_hopping(Epasc,'a','c',[0 0 0],'pz','s*','sym','E_pasc');
tb.add_hopping(-Epasc,'a','c',[1 0 0],'pz','s*','sym','E_pasc');
tb.add_hopping(-Epasc,'a','c',[0 1 0],'pz','s*','sym','E_pasc');
tb.add_hopping(Epasc,'a','c',[0 0 1],'pz','s*','sym','E_pasc');

%6th row

tb.add_hopping(0,'c','c',[0 0 0],'px','py','sym','0');
tb.add_hopping(0,'c','c',[0 0 0],'px','pz','sym','0');

tb.add_hopping(Eseapc,'c','a',[0 0 0],'px','s*','sym','E_s*apc');
tb.add_hopping(Eseapc,'c','a',[-1 0 0],'px','s*','sym','E_s*apc');
tb.add_hopping(-Eseapc,'c','a',[0 -1 0],'px','s*','sym','E_s*apc');
tb.add_hopping(-Eseapc,'c','a',[0 0 -1],'px','s*','sym','E_s*apc');

tb.add_hopping(0,'c','c',[0 0 0],'px','s*','sym','0');

%7th row
tb.add_hopping(0,'c','c',[0 0 0],'py','pz','sym','0');

tb.add_hopping(Eseapc,'c','a',[0 0 0],'py','s*','sym','E_s*apc');
tb.add_hopping(-Eseapc,'c','a',[-1 0 0],'py','s*','sym','E_s*apc');
tb.add_hopping(Eseapc,'c','a',[0 -1 0],'py','s*','sym','E_s*apc');
tb.add_hopping(-Eseapc,'c','a',[0 0 -1],'py','s*','sym','E_s*apc');

tb.add_hopping(0,'c','c',[0 0 0],'py','s*','sym','0');

%8th row
tb.add_hopping(Eseapc,'c','a',[0 0 0],'pz','s*','sym','E_s*apc');
tb.add_hopping(-Eseapc,'c','a',[-1 0 0],'pz','s*','sym','E_s*apc');
tb.add_hopping(-Eseapc,'c','a',[0 -1 0],'pz','s*','sym','E_s*apc');
tb.add_hopping(Eseapc,'c','a',[0 0 -1],'pz','s*','sym','E_s*apc');

tb.add_hopping(0,'c','c',[0 0 0],'pz','s*','sym','0');

%9th row
tb.add_hopping(0,'a','c',[0 0 0],'s*','s*','sym','0');





%g1 = 1+d1-d2-d3
%g2 = 1-d1+d2-d3
%g3 = 1-d1-d2+d3





tb.symbolic_hamiltonian()
