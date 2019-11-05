%Graphene Structure
% [A]--[B]         [A]--[B]
% (1)    \        /  (2)
%         [A]--[B] >>> this is 0
% (3)    /        \  (4)
% [A]--[B]         [A]--[B]

%From (0) to go (1) we have -a1, so we will add hopping as following
%Our translation vector R = -a1 + 0a2 so it is [-1 0]
%Hopping amplitude as t
%interaction is which atoms interact which ones, in our case 0 to 1 A atoms
%interacts with B
%add_hopping(amplitude,interaction between pairs,translation vec)
%add_hopping(t,1,2,[0 0]);  % inside 0 A and B interaction interaction
%between zero and one which A and B respectively.
%add_hopping(t,1,2,[-1 0]);  % 0 -> 1
%add_hopping(t,2,1,[1 0]);   % 0 -> 2
%add_hopping(t,2,1,[0 1]);  % 0 -> 3
%add_hopping(t,1,2,[0 -1]);  % 0 -> 4

clear; clc; close all;
%graphene consts

a0 = 1.42*1e-10;
a = a0*sqrt(3);

%primitive vectors
a1 = [a*sqrt(3)/2,  a/2,  0];
a2 = [a*sqrt(3)/2, -a/2, 0];

%a1 = [a,0,0];
%a2 = [a/2,a/2*sqrt(3),0];


tb = tightbinding('graphene',a1,a2);% Start with name and primitive vectors
tb.set_unit_cell('A',[-a0/2 0],'B',[a0/2 0]); %give unit cell atoms and their locations
%tb.set_unit_cell('A',[0, -a0/2],'B',[0, a0/2]); %give unit cell atoms and their locations


Eo = 0;
t = 0.1;
tb.add_hopping(Eo,1,1,[0 0]);  % 0 to 0
tb.add_hopping(Eo,2,2,[0 0]);  % 0 to 0

tb.add_hopping(-t,1,2,[0 0]);  % 0 to 0
tb.add_hopping(-t,2,1,[0 0]);  % 0 to 0

tb.add_hopping(-t,1,2,[-1 0]);  % 0 to 1 with -a1 vector
tb.add_hopping(-t,2,1,[1 0]);   % 0 to 2 with a1 vector
tb.add_hopping(-t,1,2,[0 -1]);  % 0 to 3
tb.add_hopping(-t,2,1,[0 1]);  % 0 to 4

%Do this, just make sure to get symmetric hamiltonian
%tb.hermitian_hamiltonian();
if(0)
lat_f = figure(3);
gp = lattice_drawer(lat_f,20,20);
atoma = gp.draw('circle blue',0,0,0.2,'Visible','off');
atomb = gp.draw('circle red',0,0,0.3,'Visible','off');
bond = gp.draw('line black',0,0,0,0,'Visible','off');
type_str.bonds = {bond};
type_str.atoms = {atoma,atomb};
tb.plot_lattice(gp,type_str);

uc_rect{1} = gp.draw('line',-2,0,0,1,'Color','red');
uc_rect{2} = gp.draw('line',0,1,2,0,'Color','red');
uc_rect{3} = gp.draw('line',2,0,0,-1,'Color','red');
uc_rect{4} = gp.draw('line',0,-1,-2,0,'Color','red');
end


if(1)

range = 2*pi/a;
precision = 100;
tb.set_kvector(-range,range,precision);
tb.calculate_band();
fig_band = figure(1);
surfaces = tb.plot_band(fig_band);
%For contour, uncomment below
%contour(tb.kvec{1},tb.kvec{2},tb.E);

b = tb.get_reciprocal_vectors();
b1 = b{1}; %*a*sqrt(3)/(2*pi);
b2 = b{2};%*a*sqrt(3)/(2*pi);
points = {};
for i = [-1,0,1]
	for j = [-1,0,1]
		points(end+1) = {i.*b1+j.*b2};
	end
end


rp = lattice_drawer(fig_band);
rp.xaxis_symmetric();
rp.yaxis_symmetric();
px = [];
py = [];
rp.draw('vector red',0,0,b1(1),b1(2));
rp.draw('vector blue',0,0,b2(1),b2(2));

for i = 1:size(points,2)
	pt = points{i};
	rp.draw('point',pt(1),pt(2),0.5);
	px(end+1) = pt(1);
	py(end+1) = pt(2);	
end

%iki vektörü sortlamk için
%[x,order] = sort(x)
%ys = y(order); y = ys; clear ys;


[vpx,vpy] = voronoi(px,py);
%voronoi(px,py);


vpx1 = vpx(1,:);
vpy1 = vpy(1,:);
[vpx1,I] = sort(vpx1);
vpy1 = vpy1(I);

good = {};

i = 1;
while( i <= size(vpx1,2) )
	pt(1) = vpx1(i);
	pt(2) = vpy1(i);
	j = i+1;
	while(j <= size(vpx1,2) )
		if(pt(1) == vpx1(j) & pt(2) == vpy1(j))
			vpx1(j) = [];
			vpy1(j) = [];
		end
		j = j+1;
	end
	i = i+1;
end
hypots = transpose(hypot(vpx1(:),vpy1(:)));

[hypots,I] = sort(hypots);
vpx1 = vpx1(I);
vpy1 = vpy1(I);



oldhyp = hypots;
hypotmin = min(hypots);
i = 1;
while (i <= size(hypots,2))
	if(abs(hypots(i)-hypotmin) > 1e-4)
		hypots(i) = [];
		vpx1(i) = [];
		vpy1(i) = [];
	else
		i = i + 1;
	end 
end

angles = atan2d(vpy1,vpx1);
[angles,I] = sort(angles);
vpx1 = vpx1(I);
vpy1 = vpy1(I);

vpx1(end+1) = vpx1(1);
vpy1(end+1) = vpy1(1);

rp.draw('line black',vpx1,vpy1,0.5*ones(1,size(vpx1,2)));

surfaces.positive_surface.Visible = 'off';
surfaces.negative_surface.Visible = 'off';

for i = 1:size(vpx1,2)
		temp = rp.draw('point yellow',vpx1(i),vpy1(i),0.5);
		t = rp.set_text(temp,sprintf("%d",i));
		t.Position(3) = 0.5;
end

view(0,90);


Gamma = [0, 0,0];
K1 = [-2*pi / (sqrt(3)*a), 2*pi / (3*a), 0];
M = [2*pi / (sqrt(3)*a), 0, 0];
K2 = [2*pi / (sqrt(3)*a), 2*pi / (3*a), 0];

k_fig = figure(2);
f2 = tb.plot_high_symmetry_points(k_fig,K1,Gamma,M,K2);
%f2 = tb.plot_high_symmetry_points(figure(2),Gamma,M,K2,Gamma);

title(f2.CurrentAxes,'FBG','Interpreter','Latex');
xlabel(f2.CurrentAxes,'$$\Gamma K M K$$','Interpreter','Latex');
xticks(f2.CurrentAxes,[0 50 100 150]);
xticklabels(f2.CurrentAxes,{'K','\Gamma','M','K'});


% band_dw = lattice_drawer(fig_band); 
% l = band_dw.draw('vector',K1(1),K1(2),1,Gamma(1),Gamma(2),1);
% %vector için
% l.MaxHeadSize = 1;
% %l.ZData = [1 1];
% l = band_dw.draw('vector',Gamma(1),Gamma(2),1,M(1),M(2),1);
% %l.ZData = [1 1];
% l.MaxHeadSize = 1;
% l = band_dw.draw('vector',M(1),M(2),1,K2(1),K2(2),1);
% %l.ZData = [1 1];
% l.MaxHeadSize = 1;
% l = band_dw.draw('vector',K2(1),K2(2),1,K1(1),K1(2),1);
% %l.ZData = [1 1];
% l.MaxHeadSize = 1;

end


%rp.draw('vector black',0,0,b1(1)+b2(1),b1(2)+b2(2));
%rp.draw('vector yellow',0,0,b1(1)-b2(1),b1(2)-b2(2));






%gp_byhand;
%f.Color = [0.5 0.5 0.2];

%K1 = [-4*pi / (3*sqrt(3)*a), 0,0];
%Gamma = [0, 0, 0];
%gp = lattice_drawer(f);
%ll = gp.draw('vector',K1(1),K1(2),Gamma(1),Gamma(2));
%ll.MaxHeadSize = 1;
%ll.LineWidth = 2;
%M = [0, 2*pi / (3*a),0];
%ll = gp.draw('vector',Gamma(1),Gamma(2),M(1),M(2));
%ll.MaxHeadSize = 1;
%ll.LineWidth = 2;
%K2 = [2*pi / (3*sqrt(3)*a), 2*pi / (3*a),0];
%ll = gp.draw('vector',M(1),M(2),K2(1),K2(2));
%ll.LineWidth = 2;
