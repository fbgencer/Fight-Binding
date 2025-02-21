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

%Soc params
soa=.3787/3;
soc=.0129/3;

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
tb.spin_orbit_coupling('true');

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


%%first row
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


%Spin-Orbit part
tb.add_soc(-i*soa,'a+','a+','px','py','sym',-i*sym('Ga','real'));
tb.add_soc(soa,'a+','a-','px','pz','sym',sym('Ga','real'));
tb.add_soc(-i*soa,'a+','a-','py','pz','sym',-i*sym('Ga','real'));
tb.add_soc(-soa,'a+','a-','pz','px','sym',-sym('Ga','real'));
tb.add_soc(i*soa,'a+','a-','pz','py','sym',i*sym('Ga','real'));
tb.add_soc(i*soa,'a-','a-','px','py','sym',i*sym('Ga','real'));

tb.add_soc(-i*soc,'c+','c+','px','py','sym',-i*sym('Gc','real'));
tb.add_soc(soc,'c+','c-','px','pz','sym',sym('Gc','real'));
tb.add_soc(-i*soc,'c+','c-','py','pz','sym',i*sym('Gc','real'));
tb.add_soc(-soc,'c+','c-','pz','px','sym',-sym('Gc','real'));
tb.add_soc(i*soc,'c+','c-','pz','py','sym',i*sym('Gc','real'));
tb.add_soc(i*soc,'c-','c-','px','py','sym',i*sym('Gc','real'));

%S = tb.symbolic_hamiltonian();
%return

if(0)
lat_f = figure("Name","Lattice Figure");
gp = lattice_drawer(lat_f,3,3,3);


atoma = gp.draw('sphere rgb:33FFFF',0,0,0,0.1,'Visible','off');
atomb = gp.draw('point black',0,0,0.1,'Visible','off');
bond = gp.draw('line rgb:FF55FF',0,0,0,0,'Visible','off','LineWidth',0.01,'LineStyle','--');
bonds = {bond};
atoms = {atoma,atomb};
from_to = 1:2;
tb.plot_lattice(gp,"bonds",bonds,"atoms",atoms,"x",from_to,"y",from_to,"z",from_to);
gp.set_xlabel('X');
gp.set_ylabel('Y');
gp.set_zlabel('Z');
grid;

an = a*1e10;
gp.draw('cuboid',an,an,an,an,an,an,'FaceColor','None','LineWidth',1);
%a1 = a1.*1e10;
%a2 = a2.*1e10;
%a3 = a3.*1e10;
%v = gp.draw('vector black',0.5,0.5,0.5,a1(1),a1(2),a1(3),'MaxHeadSize',0.6);
%v = gp.draw('vector black',0.5,0.5,0.5,a2(1),a2(2),a2(3),'MaxHeadSize',0.6);
%v = gp.draw('vector black',0.5,0.5,0.5,a3(1),a3(2),a3(3),'MaxHeadSize',0.6);



end

if(0)
range = 2*pi/a; 
precision = 15;
k = tb.set_kvector(-range,range,precision);
fig_band = figure("Name","Energy Band Figure");
surfaces = tb.plot_energy_band(fig_band,k,'surface','EdgeColor','None');
end

if(1)
%both points are working
G = [0 0 0];
X = [2*pi/a 0 0];
L = [pi/a pi/a pi/a];
precision =  100;

X = {2*pi/a*[1 0 0],'X'};
G = {[0 0 0],'$$\Gamma$$'};
L = {(pi/a)*[1 1 1],'L'};

fig_hsym = figure("Name","High Symmetry Points Figure");
tb.plot_high_symmetry_points(fig_hsym,precision,L,G,X);

%fig_hsym.CurrentAxes.YLim = [-3 3];

end

if(0)
Gamma = [0 0 0];
X = [2*pi/a 0 0];
L = [pi/a pi/a pi/a];


fig_hsym = figure("Name","High Symmetry Points Figure");
plts = tb.plot_high_symmetry_points(fig_hsym,500,L,Gamma,X);
xlabel(fig_hsym.CurrentAxes,'$$L \Gamma X$$','Interpreter','Latex');
xticks(fig_hsym.CurrentAxes,[0 100 200]);
xticklabels(fig_hsym.CurrentAxes,{'L','\Gamma','X'});
ylabel(fig_hsym.CurrentAxes,'$$Energy(eV)$$','Interpreter','Latex')
grid();
%fig_hsym.CurrentAxes.YLim = [-3 3];
%for ix = 1:numel(plts), plts{ix}.Color = 'red'; end
end


%g1 = 1+d1-d2-d3
%g2 = 1-d1+d2-d3
%g3 = 1-d1-d2+d3

toc;

if(0)
fig = figure("Name","Datta");

soa=.3787/3;soc=.0129/3;
Esa=-8.3431;Epa=1.0414;Esc=-2.6569;Epc=3.6686;Esea=8.5914;
Esec=6.7386;
Vss=-6.4513;Vxx=1.9546;Vxy=5.0779;Vsapc=4.4800;Vpasc=5.7839;Vseapc=4.8422;Vpasec=4.8077;
d1=[1 1 1]/4;d2=[1 -1 -1]/4;d3=[-1 1 -1]/4;d4=[-1 -1 1]/4;
d1=[0 0 0]/2;d2=[0 -1 -1]/2;d3=[-1 0 -1]/2;d4=[-1 -1 0]/2;

for iter = 1:2
if iter == 1
l=1;m=1;n=1;kmax=pi;Nt=101;%L-direction
else
l=1;m=0;n=0;kmax=2*pi;Nt=101;%X-direction
end


for Nk=1:Nt
k=[l m n]*kmax*(Nk-1)/(Nt-1);
	p1=exp(i*sum(k.*d1));p2=exp(i*sum(k.*d2));
	p3=exp(i*sum(k.*d3));p4=exp(i*sum(k.*d4));
		g0=(p1+p2+p3+p4)/4;g1=(p1+p2-p3-p4)/4;
		g2=(p1-p2+p3-p4)/4;g3=(p1-p2-p3+p4)/4;

h=[Esa/2 Vss*g0 0 0 0 Vsapc*g1 Vsapc*g2 Vsapc*g3 0 0;
      0	Esc/2	-Vpasc*conj(g1)	-Vpasc*conj(g2)      -Vpasc*conj(g3)    0 0 0 0 0;
      0	0	Epa/2	0	0	Vxx*g0 Vxy*g3	Vxy*g2  0  -Vpasec*g1;
      0	0	0	Epa/2	0	Vxy*g3 Vxx*g0	Vxy*g1  0  -Vpasec*g2;
      0	0	0	0	Epa/2	Vxy*g2 Vxy*g1	Vxx*g0  0  -Vpasec*g3;
      0	0	0	0	0	Epc/2	0	0	Vseapc*(g1)	0;
      0	0	0	0	0	0	Epc/2	0	Vseapc*(g2)	0;
      0	0	0	0	0	0	0	Epc/2	Vseapc*(g3)	0;
      0	0	0	0	0	0	0	0	Esea/2	0;
      0	0	0	0	0	0	0	0	0	Esec/2];

H=[h+h'	zeros(10);
     zeros(10)	h+h'];

hso=zeros(20);
hso(3,4)=-i*soa;hso(3,15)=soa;
hso(4,15)=-i*soa;
hso(5,13)=-soa;hso(5,14)=i*soa;
hso(6,7)=-i*soc;hso(6,18)=soc;
hso(7,18)=-i*soc;
hso(8,16)=-soc;hso(8,17)=i*soc;
hso(13,14)=i*soa;
hso(16,17)=i*soc;
Hso=hso+hso';

eiglst = eig(H+Hso);
        E(Nk,:) = sort(real(eiglst));
		X(Nk)=-(Nk-1)/(Nt-1);%L-direction
		X1(Nk)=(Nk-1)/(Nt-1);%X-direction
end

hold on
if iter == 1
h=plot(X,E,'b');
else
h=plot(X1,E,'b');

end
axis([-1 1 -3 3])
xlabel('k (as fraction of maximum value)--->')
ylabel('Energy (eV) ---> ')
grid on
end

end

