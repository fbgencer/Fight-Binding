% NbSe2
clear all; clc; close all;
a = 1;

%primitive vectors
a1 = [3.488132 0.000000 0.000000];
a2 = [-1.744066 3.020811 0.000000];
a3 = [0.000000 0.000000 13.873137];



tb = tightbinding(3,a1,a2,a3);% Start with name and primitive vectors
Nb1 =[-0.000002    2.013875    3.468284];
Nb2=[1.744068    1.006936    10.404853];
Se1=[-0.000002    2.013875    12.095073];
Se2=[1.744068    1.006936    1.778065];
Se3=[1.744068    1.006936    5.158504];
Se4=[-0.000002    2.013875    8.714634];

tb.set_unit_cell('Nb',Nb1,'Nb',Nb2,...
	'Se',Se1,'Se',Se2,'Se',Se3,'Se',Se4);

tb.set_orbital('d1,d2,d3,d4');
tb.set_metric_unit('A');
tb.spin_orbit_coupling('false');
tb.add_hrfile('NbSe2_hr_7007.dat');
tb.set_fermi_level(9);



if(0)
% 0.0000 0.0000 0.0000 20 !Gamma
% 0.5000 0.0000 0.0000 20 !M
% 0.3300 0.3300 0.0000 20 !K
% 0.0000 0.0000 0.0000 20 !Gamma
% 0.5000 0.0000 0.5000 20 !L
G = {[0 0 0],'$$\Gamma$$'};
M = {2*pi*[0.5 0 0],'M'};
K = {2*pi*[0.33 0.33 0],'K'};
L = {2*pi*[0.5 0 0.5],'L'};

precision = [100,30,100];

fig_hsym = figure("Name","High Symmetry Points Figure");
tb.plot_high_symmetry_points(fig_hsym,precision,G,M,K,G);
ax = fig_hsym.CurrentAxes;
ax.YLim = [-6 5];
%ax.YTick = -5:1:5;
%ax.YTickLabel = {'-5','-4','-3','-2','-1','0','1','2','3','4','5'}
%,W,L,K,G
end




if(1)
lat_f = figure("Name","Lattice Figure");
gp = lattice_drawer(lat_f,20,20,30);
gp.xaxis_symmetric();
gp.yaxis_symmetric();
gp.zaxis_symmetric();


%Gallium covalent radius = 1.26 ang, 1.36 atomic
%As covalent radius = 1.19 ang, 1.14 atomic

atoma = gp.draw('sphere rgb:4589e4',0,0,0,1.36/2,'Visible','off');
atomb = gp.draw('sphere black',0,0,0,1.14/2,'Visible','off');
bond = gp.draw('line red',0,0,0,0,'Visible','off','LineWidth',4,'LineStyle','-');
%rgb:FF55FF
bonds = {bond};
atoms = {atoma,atomb};
from_to = 10:40;
tb.plot_lattice(gp,"bonds",bonds,"atoms",atoms);
%,"z",-6:6,"x",-3:3,"y",-3:3
%legend({'Nb','Se'},'Location','southwest') 



orj = [0 0 0];

v1 = orj;
v2 = a1;
px = [v1(1) v2(1)];
py = [v1(2) v2(2)];
pz = [v1(3) v2(3)];
gp.draw('line green',px,py,pz,'LineWidth',3);


v1 = orj;
v2 = a2;
px = [v1(1) v2(1)];
py = [v1(2) v2(2)];
pz = [v1(3) v2(3)];
gp.draw('line green',px,py,pz,'LineWidth',3);


v1 = a2;
v2 = a1+a2;
px = [v1(1) v2(1)];
py = [v1(2) v2(2)];
pz = [v1(3) v2(3)];
gp.draw('line green',px,py,pz,'LineWidth',3);



v1 = a1+a2;
v2 = a1;
px = [v1(1) v2(1)];
py = [v1(2) v2(2)];
pz = [v1(3) v2(3)];
gp.draw('line green',px,py,pz,'LineWidth',3);

%%%%%%%%%%%%%%%%%%%%%%
s = orj;
x = [];
y = [];
z = [];

for k = 0:1
	for i = 0:1
		for j = 0:1
			v1 = i*a1+j*a2+k*a3
			%v2 = s;
			px = [v1(1) v2(1)];
			py = [v1(2) v2(2)];
			pz = [v1(3) v2(3)];
			gp.draw('point',v1(1),v1(2),v1(3));
			x(end+1) = v1(1);
			y(end+1) = v1(2);
			z(end+1) = v1(3);
			%gp.draw('line green',px,py,pz,'LineWidth',3);
			%s = v1;
		end
	end
end


%an = a;
%gp.draw('cuboid',an,an,an,an,an,an,'FaceColor','None','LineWidth',1);

end