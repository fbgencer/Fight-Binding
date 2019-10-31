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

tb = tightbinding('graphene',a1,a2);% Start with name and primitive vectors
tb.set_unit_cell('A',[-a0/2 0],'B',[a0/2 0]); %give unit cell atoms and their locations

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
atoma = gp.draw('circle blue',0,0,0.3,'Visible','off');
atomb = gp.draw('circle red',0,0,0.3,'Visible','off');
bond = gp.draw('line black',0,0,0,0,'Visible','off');
type_str.bond = {bond};
type_str.atom = {atoma,atomb};
tb.plot_lattice(gp,type_str);

uc_rect{1} = gp.draw('line',-2,0,0,1,'Color','red');
uc_rect{2} = gp.draw('line',0,1,2,0,'Color','red');
uc_rect{3} = gp.draw('line',2,0,0,-1,'Color','red');
uc_rect{4} = gp.draw('line',0,-1,-2,0,'Color','red');
end



range = 5*pi/a;
precision = 100;
tb.set_kvector(-range,range,precision);
tb.calculate_band();
fig_band = figure(1);
%tb.plot_band(fig_band);
contour(tb.kvec{1},tb.kvec{2},tb.E);


Gamma = [0, 0,0];
K1 = [-2*pi / (sqrt(3)*a), 2*pi / (3*a), 0];
M = [2*pi / (sqrt(3)*a), 0, 0];
K2 = [2*pi / (sqrt(3)*a), 2*pi / (3*a), 0];

%k_fig = figure(2);
%f2 = tb.plot_high_symmetry_points(k_fig,K1,Gamma,M,K2);
%f2 = tb.plot_high_symmetry_points(figure(2),Gamma,M,K2,Gamma);

%title(f2.CurrentAxes,'FBG','Interpreter','Latex');
%xlabel(f2.CurrentAxes,'$$\Gamma K M K$$','Interpreter','Latex');
%xticks(f2.CurrentAxes,[0 50 100 150]);
%xticklabels(f2.CurrentAxes,{'K','\Gamma','M','K'});


band_dw = lattice_drawer(fig_band); 
l = band_dw.draw('vector',K1(1),K1(2),Gamma(1),Gamma(2));
%vector için
l.MaxHeadSize = 1;
%l.ZData = [1 1];
l = band_dw.draw('vector',Gamma(1),Gamma(2),M(1),M(2));
%l.ZData = [1 1];
l.MaxHeadSize = 1;
l = band_dw.draw('vector',M(1),M(2),K2(1),K2(2));
%l.ZData = [1 1];
l.MaxHeadSize = 1;
l = band_dw.draw('vector',K2(1),K2(2),K1(1),K1(2));
%l.ZData = [1 1];
l.MaxHeadSize = 1;



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
