% Gaas

clear; clc; close all;
tic;


Esa=-8.3431;
Epa=1.0414;
Esc=-2.6569;
Epc=3.6686;
Esea=8.5914;
Esec=6.7386;
Vss=-6.4513;
Vxx=1.9546;
Vxy=5.0779;
Vsapc=4.4800;
Vpasc=5.7839;
Vseapc=4.8422;
Vpasec=4.8077;

a = 1e-10;

%primitive vectors
a1 = (a/2).*[0, 1, 1];
a2 = (a/2).*[1, 0, 1];
a3 = (a/2).*[1, 1, 0];

g0 = [1; 1; 1; 1]/4;
g1 = [1; 1;-1;-1]/4;
g2 = [1;-1; 1;-1]/4;
g3 = [1;-1;-1; 1]/4;

d0 = [0 0 0]; % phase terms
d1 = [1 0 0];
d2 = [0 1 0];
d3 = [0 0 1];

d = [d0;d1;d2;d3];

tb = tightbinding(3,a1,a2,a3);% Start with name and primitive vectors

tb.set_unit_cell('a',[0 0 0],'c',[1 1 1].*(a/4) );
tb.set_orbital('3s','px,py,pz','4s');

%Diagonal terms
tb.add_hopping(Esa ,'a','a',d0,'3s','3s','sym','E_sa');
tb.add_hopping(Esc ,'c','c',d0,'3s','3s','sym','E_sc');
tb.add_hopping(Epa ,'a','a',d0,'px','px','sym','E_pa');
tb.add_hopping(Epa ,'a','a',d0,'py','py','sym','E_pa');
tb.add_hopping(Epa ,'a','a',d0,'pz','pz','sym','E_pa');
tb.add_hopping(Epc ,'c','c',d0,'px','px','sym','E_pc');
tb.add_hopping(Epc ,'c','c',d0,'py','py','sym','E_pc');
tb.add_hopping(Epc ,'c','c',d0,'pz','pz','sym','E_pc');
tb.add_hopping(Esea,'a','a',d0,'4s','4s','sym','E_s*a');
tb.add_hopping(Esec,'c','c',d0,'4s','4s','sym','E_s*c');


% %%first row
%Vss*g0 0 0 0 Vsapc*g1 Vsapc*g2 Vsapc*g3 0 0;
tb.add_hopping(Vss.*g0		,'a','c',d ,'3s','3s','sym','E_ss');
tb.add_hopping(0			,'a','a',d0,'3s','px','sym','0');
tb.add_hopping(0			,'a','a',d0,'3s','py','sym','0');
tb.add_hopping(0			,'a','a',d0,'3s','pz','sym','0');
tb.add_hopping(Vsapc.*(g1)	,'a','c',d ,'3s','px','sym','E_sapc');
tb.add_hopping(Vsapc.*(g2)	,'a','c',d ,'3s','py','sym','E_sapc');
tb.add_hopping(Vsapc.*(g3)	,'a','c',d ,'3s','pz','sym','E_sapc');
tb.add_hopping(0			,'a','a',d0,'3s','4s','sym','0');
tb.add_hopping(0			,'a','c',d0,'3s','4s','sym','0');

%2nd row
%-Vpasc*conj(g1)	-Vpasc*conj(g2)	-Vpasc*conj(g3) 0 0 0 0 0;
tb.add_hopping(-Vpasc.*g1,'c','a',-d,'3s','px','sym','E_pasc');
tb.add_hopping(-Vpasc.*g2,'c','a',-d,'3s','py','sym','E_pasc');
tb.add_hopping(-Vpasc.*g3,'c','a',-d,'3s','pz','sym','E_pasc');
tb.add_hopping(0		 ,'c','c',d0,'3s','px','sym','0');
tb.add_hopping(0		 ,'c','c',d0,'3s','py','sym','0');
tb.add_hopping(0		 ,'c','c',d0,'3s','pz','sym','0');
tb.add_hopping(0		 ,'c','a',d0,'3s','4s','sym','0');
tb.add_hopping(0		 ,'c','c',d0,'3s','4s','sym','0');

%3rd row
%0	0	Vxx*g0 Vxy*g3	Vxy*g2  0  -Vpasec*g1;
tb.add_hopping(0			,'a','a',d0,'px','py','sym','0');
tb.add_hopping(0			,'a','a',d0,'px','pz','sym','0');
tb.add_hopping(Vxx.*g0		,'a','c',d ,'px','px','sym','E_xx');
tb.add_hopping(Vxy.*g3		,'a','c',d ,'px','py','sym','E_xy');
tb.add_hopping(Vxy.*g2		,'a','c',d ,'px','pz','sym','E_xy');
tb.add_hopping(0			,'a','a',d0,'px','4s','sym','0');
tb.add_hopping(-Vpasec.*g1	,'a','c',d ,'px','4s','sym','E_pasc');

%4th row
%0	Vxy*g3 Vxx*g0	Vxy*g1  0  -Vpasec*g2;
tb.add_hopping(0			,'a','a',d0,'py','pz','sym','0');
tb.add_hopping(Vxy.*g3		,'a','c',d ,'py','px','sym','E_xy');
tb.add_hopping(Vxx.*g0		,'a','c',d ,'py','py','sym','E_xx');
tb.add_hopping(Vxy.*g1		,'a','c',d ,'py','pz','sym','E_xy');
tb.add_hopping(0			,'a','a',d0,'py','4s','sym','0');
tb.add_hopping(-Vpasec.*g2	,'a','c',d ,'py','4s','sym','E_pasc');

%5th row
%Vxy*g2 Vxy*g1	Vxx*g0  0  -Vpasec*g3;
tb.add_hopping(Vxy.*g2		,'a','c',d ,'pz','px','sym','E_xy');
tb.add_hopping(Vxy.*g1		,'a','c',d ,'pz','py','sym','E_xy');
tb.add_hopping(Vxx.*g0		,'a','c',d ,'pz','pz','sym','E_xx');
tb.add_hopping(0			,'a','a',d0,'pz','4s','sym','0');
tb.add_hopping(-Vpasec.*g3	,'a','c',d ,'pz','4s','sym','E_pasc');

% %6th row
%0	0	Vseapc*(g1)	0;
tb.add_hopping(0			,'c','c',d0,'px','py','sym','0');
tb.add_hopping(0			,'c','c',d0,'px','pz','sym','0');
tb.add_hopping(Vseapc.*g1	,'c','a',d ,'px','4s','sym','E_s*apc');
tb.add_hopping(0			,'c','c',d0,'px','4s','sym','0');

%7th row
%0	Vseapc*(g2)	0;
tb.add_hopping(0			,'c','c',d0,'py','pz','sym','0');
tb.add_hopping(Vseapc.*g2	,'c','a',d ,'py','4s','sym','E_s*apc');
tb.add_hopping(0			,'c','c',d0,'py','4s','sym','0');

%8th row
%Vseapc*(g3)	0;
tb.add_hopping(Vseapc.*g3,'c','a',d ,'pz','4s','sym','E_s*apc');
tb.add_hopping(0		 ,'c','c',d0,'pz','4s','sym','0');

%9th row
%0
tb.add_hopping(0,'a','c',d0,'4s','4s','sym','0');


%sym_ham = tb.symbolic_hamiltonian();

if(1)
lat_f = figure("Name","Lattice Figure");
gp = lattice_drawer(lat_f,3,3,3);


atoma = gp.draw('point blue',0,0,0.25,'Visible','off','LineWidth',5);
atomb = gp.draw('point black',0,0,0.1,'Visible','off');
bond = gp.draw('line rgb:FF33FF',0,0,0,0,'Visible','off','LineWidth',0.05);
bonds = {bond};
atoms = {atoma,atomb};
from_to = 1:2;
tb.plot_lattice(gp,"bonds",bonds,"atoms",atoms,"x",from_to,"y",from_to,"z",from_to);
gp.set_xlabel('X');
gp.set_ylabel('Y');
gp.set_zlabel('Z');
grid;

a = a*1e10;
gp.draw('cuboid',a,a,a,a,a,a,'FaceColor','None','LineWidth',1);

%gp.xaxis_symmetric

end

if(0)
range = 2*pi/a; 
precision = 15;
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
grid();
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
	p1 = 1;
	p2=exp(-i*dot(k,d1));
	p3=exp(-i*dot(k,d2));
	p4=exp(-i*dot(k,d3));
	g0=(p1+p2+p3+p4)/4;
	g1=(p1+p2-p3-p4)/4;
	g2=(p1-p2+p3-p4)/4;
	g3=(p1-p2-p3+p4)/4;

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
eiglst = eig(H);
        %eiglst = sum(D);
        E(Nk,:) = ((eiglst));
        
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

