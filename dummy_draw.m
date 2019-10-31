close all; clear; clc;
gp = lattice_drawer(figure(1),30,30);

%o1 = gp.draw('line',3,3,3,10);
%o1 = gp.draw('line',3,10,25,20);
%o1 = gp.draw('line',25,20,3,3);


gp.set_title('$Graphene\hspace{1mm}C_{6}$','Interpreter','latex');
gp.set_xlabel('$X-Coordinates$','Interpreter','latex');
gp.set_ylabel('$Y-Coordinates$','Interpreter','latex');

%dumy = gp.draw('rect',5,5,3,4,'FaceColor','blue','Visible','off');
%a = gp.copy_to(dumy,10,10,10,2,'Visible','on');

cr = gp.draw('circle black',10,10,3);
copy_cr =  gp.copy_to(cr,20,10,'EdgeColor','black','FaceColor','red');

rc = gp.draw('rect blue',13.5,18.5,3,3);
l = gp.draw('line',15,0,15,30,'Color','red');
bond = gp.draw('line black',0,0,0,0,'Visible','off');
lc = gp.copy_to(bond,0,20,30,20,'Visible','on');

%gp.set_text(l,'FBGENCER');