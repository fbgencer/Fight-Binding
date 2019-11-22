% Gaas

clear; clc; close all;
tic;

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
Epasec=4.8077;


%Esa=-8.3431;Epa=1.0414;Esc=-2.6569;Epc=3.6686;Esea=8.5914;Esec=6.7386;
%Vss=-6.4513;Vxx=1.9546;Vxy=5.0779;Vsapc=4.4800;Vpasc=5.7839;Vseapc=4.8422;
%Vpasec=4.8077;

a = 1;

%primitive vectors
a1 = (a/2).*[0, 1, 1];
a2 = (a/2).*[1, 0, 1];
a3 = (a/2).*[1, 1, 0];

g0 = [1; 1; 1; 1];
g1 = [1; 1;-1;-1];
g2 = [1;-1; 1;-1];
g3 = [1;-1;-1; 1];

d0 = [0 0 0]; % phase terms
d1 = [1 0 0];
d2 = [0 1 0];
d3 = [0 0 1];

d = [d0;d1;d2;d3];

tb = tightbinding(3,a1,a2,a3);% Start with name and primitive vectors

tb.set_unit_cell('a',[0 0 0],'c',[1 1 1].*a/4);
tb.set_orbital('3s','px,py,pz','4s','d1,d2,d3');

%Diagonal terms
%tb.add_hopping(Esa ,'a','a',d0,'3s','px','sym','E_sa');
tb.add_hopping(Esc ,'a','c',d0,'3s','d3','sym','E_sc');
tb.symbolic_hamiltonian()
% tb.add_hopping(Epa ,'a','a',d0,'3px','3px','sym','E_pa');
% tb.add_hopping(Epa ,'a','a',d0,'3py','3py','sym','E_pa');
% tb.add_hopping(Epa ,'a','a',d0,'3pz','3pz','sym','E_pa');
% tb.add_hopping(Epc ,'c','c',d0,'3px','3px','sym','E_pc');
% tb.add_hopping(Epc ,'c','c',d0,'3py','3py','sym','E_pc');
% tb.add_hopping(Epc ,'c','c',d0,'3pz','3pz','sym','E_pc');
% tb.add_hopping(Esea,'a','a',d0,'4s','4s','sym','E_s*a');
% tb.add_hopping(Esec,'c','c',d0,'4s','4s','sym','E_s*c');


% %%first row
% %Vss*g0 0 0 0 Vsapc*g1 Vsapc*g2 Vsapc*g3 0 0;
% tb.add_hopping(Ess.*g0,'a','c',d,'s','s','sym','E_ss');

% tb.add_hopping(0,'a','a',d0,'s','px','sym','0');
% tb.add_hopping(0,'a','a',d0,'s','py','sym','0');
% tb.add_hopping(0,'a','a',d0,'s','pz','sym','0');

% tb.add_hopping(Esapc.*(g1),'a','c',d,'s','px','sym','E_sapc');
% tb.add_hopping(Esapc.*(g2),'a','c',d,'s','px','sym','E_sapc');
% tb.add_hopping(Esapc.*(g3),'a','c',d,'s','px','sym','E_sapc');

% tb.add_hopping(0,'a','a',d0,'s','s*','sym','0');
% tb.add_hopping(0,'a','c',d0,'s','s*','sym','0');

% %2nd row
% %-Vpasc*conj(g1)	-Vpasc*conj(g2)	-Vpasc*conj(g3) 0 0 0 0 0;
% tb.add_hopping(-Epasc.*(-g1),'c','a',d,'s','px','sym','E_pasc');
% tb.add_hopping(-Epasc.*(-g2),'c','a',d,'s','py','sym','E_pasc');
% tb.add_hopping(-Epasc.*(-g3),'c','a',d,'s','pz','sym','E_pasc');

% tb.add_hopping(0,'c','c',d0,'s','px','sym','0');
% tb.add_hopping(0,'c','c',d0,'s','py','sym','0');
% tb.add_hopping(0,'c','c',d0,'s','pz','sym','0');

% tb.add_hopping(0,'c','a',d0,'s','s*','sym','0');
% tb.add_hopping(0,'c','c',d0,'s','s*','sym','0');

% %3rd row
% %0	0	Vxx*g0 Vxy*g3	Vxy*g2  0  -Vpasec*g1;
% tb.add_hopping(0,'a','a',d0,'px','py','sym','0');
% tb.add_hopping(0,'a','a',d0,'px','pz','sym','0');

% tb.add_hopping(Exx.*g0,'a','c',d,'px','px','sym','E_xx');
% tb.add_hopping(Exy.*g3,'a','c',d,'px','py','sym','E_xy');
% tb.add_hopping(Exy.*g2,'a','c',d,'px','pz','sym','E_xy');

% tb.add_hopping(0,'a','a',d0,'px','s*','sym','0');

% %-Epasec
% tb.add_hopping(-Epasec.*g1,'a','c',d,'px','s*','sym','E_pasc');

% %4th row
% %0	Vxy*g3 Vxx*g0	Vxy*g1  0  -Vpasec*g2;
% tb.add_hopping(0,'a','a',d0,'py','pz','sym','0');

% tb.add_hopping(Exy.*g3,'a','c',d,'py','px','sym','E_xy');
% tb.add_hopping(Exx.*g0,'a','c',d,'py','py','sym','E_xx');
% tb.add_hopping(Exy.*g1,'a','c',d,'py','pz','sym','E_xy');

% tb.add_hopping(0,'a','a',d0,'py','s*','sym','0');

% %-Epasec
% tb.add_hopping(-Epasec.*g2,'a','c',d,'py','s*','sym','E_pasc');

% %5th row
% %Vxy*g2 Vxy*g1	Vxx*g0  0  -Vpasec*g3;
% tb.add_hopping(Exy.*g2,'a','c',d,'pz','px','sym','E_xy');
% tb.add_hopping(Exy.*g1,'a','c',d,'pz','py','sym','E_xy');
% tb.add_hopping(Exx.*g0,'a','c',d,'pz','pz','sym','E_xx');

% tb.add_hopping(0,'a','a',d0,'pz','s*','sym','0');

% %-Epasec
% tb.add_hopping(-Epasec.*g3,'a','c',d,'pz','s*','sym','E_pasc');

% % %6th row
% %0	0	Vseapc*(g1)	0;
% tb.add_hopping(0,'c','c',d0,'px','py','sym','0');
% tb.add_hopping(0,'c','c',d0,'px','pz','sym','0');

% tb.add_hopping(Eseapc.*g1,'c','a',d,'px','s*','sym','E_s*apc');

% tb.add_hopping(0,'c','c',d0,'px','s*','sym','0');

% %7th row
% %0	Vseapc*(g2)	0;
% tb.add_hopping(0,'c','c',d0,'py','pz','sym','0');
% tb.add_hopping(Eseapc.*g2,'c','a',d,'px','s*','sym','E_s*apc');
% tb.add_hopping(0,'c','c',d0,'py','s*','sym','0');

% %8th row
% %Vseapc*(g3)	0;
% tb.add_hopping(Eseapc.*g3,'c','a',d,'px','s*','sym','E_s*apc');
% tb.add_hopping(0,'c','c',d0,'pz','s*','sym','0');

% %9th row
% %0
% tb.add_hopping(0,'a','c',d0,'s*','s*','sym','0');


%sym_ham = tb.symbolic_hamiltonian();

if(0)
lat_f = figure("Name","Lattice Figure");
gp = lattice_drawer(lat_f,6,6,6);
atoma = gp.draw('point red',0,0,0.25,'Visible','off');
atomb = gp.draw('point black',0,0,0.25,'Visible','off');
bond = gp.draw('line black',0,0,0,0,'Visible','off');
bonds = {bond};
atoms = {atoma,atomb};
tb.plot_lattice(gp,"bonds",bonds,"atoms",atoms,"x",-3:3,"y",-1:1,"z",-1:1);
end

if(0)
range = 2*pi/a; 
precision = 10;
k = tb.set_kvector(-range,range,precision);
fig_band = figure("Name","Energy Band Figure");
surfaces = tb.plot_energy_band(fig_band,k,'surface','EdgeColor','None');
end


if(0)
Gamma = [0 0 0];
X = [2*pi/a 0 0];
L = [pi/a pi/a pi/a];


fig_hsym = figure("Name","High Symmetry Points Figure");
tb.plot_high_symmetry_points(fig_hsym,L,Gamma,X);
xlabel(fig_hsym.CurrentAxes,'$$L \Gamma X$$','Interpreter','Latex');
xticks(fig_hsym.CurrentAxes,[0 100 200]);
xticklabels(fig_hsym.CurrentAxes,{'L','\Gamma','X'});

end


%g1 = 1+d1-d2-d3
%g2 = 1-d1+d2-d3
%g3 = 1-d1-d2+d3

toc;

if(0)
fig = figure("Name","Datta");

Esa=-8.3431;Epa=1.0414;Esc=-2.6569;Epc=3.6686;Esea=8.5914;Esec=6.7386;
Vss=-6.4513;Vxx=1.9546;Vxy=5.0779;Vsapc=4.4800;Vpasc=5.7839;Vseapc=4.8422;
Vpasec=4.8077;

%Either of the following choices for d1,d2,d3 and d4 should give the same result.
%d1=[1 1 1]/4;d2=[1 -1 -1]/4;d3=[-1 1 -1]/4;d4=[-1 -1 1]/4;
%d1=[0 0 0]/2;d2=[0 -1 -1]/2;d3=[-1 0 -1]/2;d4=[-1 -1 0]/2;
d0=[0 0 0]/2;d1=[0 1 1]/2;d2=[1 0 1]/2;d3=[1 1 0]/2;


for iter = 1:2

if(iter == 1)
l=1;m=1;n=1;kmax=pi;Nt=100;%L-direction
else
l=1;m=0;n=0;kmax=2*pi;Nt=100;%X-direction
end

for Nk=1:Nt
k=[l m n]*kmax*(Nk-1)/(Nt-1);
	p1 = 1; p2=exp(-i*sum(k.*d1));
	p3=exp(-i*sum(k.*d2));p4=exp(-i*sum(k.*d3));
		g0=(p1+p2+p3+p4)/4;g1=(p1+p2-p3-p4)/4;
		g2=(p1-p2+p3-p4)/4;g3=(p1-p2-p3+p4)/4;

h=[	Esa/2 	Vss*g0 	0 		0 0 Vsapc*g1 Vsapc*g2 Vsapc*g3 0 0;
    0		Esc/2	-Vpasc*conj(g1)	-Vpasc*conj(g2)	-Vpasc*conj(g3)     0 0 0 0 0;
    0		0		Epa/2	0	0	Vxx*g0 Vxy*g3	Vxy*g2  0  -Vpasec*g1;
    0		0		0		Epa/2	0	Vxy*g3 Vxx*g0	Vxy*g1  0  -Vpasec*g2;
    0		0		0			0	Epa/2	Vxy*g2 Vxy*g1	Vxx*g0  0  -Vpasec*g3;
    0		0		0				0	0	Epc/2	0	0	Vseapc*(g1)	0;
    0		0		0					0	0	0	Epc/2	0	Vseapc*(g2)	0;
    0		0		0						0	0	0	0	Epc/2	Vseapc*(g3)	0;
    0		0		0							0	0	0	0	0	Esea/2	0;
    0		0		0								0	0	0	0	0	0	Esec/2];
H=h+h';
[V,D]=eig(H);
        eiglst = sum(D);
        E(Nk,:) = (real(eiglst));
        
        if(iter == 1)
			X(Nk)=-(Nk-1)/(Nt-1);%L-direction
		else
			X1(Nk)=(Nk-1)/(Nt-1);%X-direction
		end
end

hold on
if(iter == 1)
	plt = plot(X,E);
else
	plt = plot(X1,E);
end
set(plt,'linewidth',[1])
set(gca,'Fontsize',[12])
xlabel('k (as fraction of maximum value)--->')
ylabel('Energy (eV) ---> ')
grid on
%Note: X-axis from 0 to +1 represents the -X direction
%while the section from 0 to –1 represents the -L direction
end

end


