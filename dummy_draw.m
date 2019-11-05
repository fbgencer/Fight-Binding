close all; clear; clc;
gp = lattice_drawer(figure(1),30,30);

%o1 = gp.draw('line',3,3,3,10);
%o1 = gp.draw('line',3,10,25,20);
%o1 = gp.draw('line',25,20,3,3);


gp.set_title("$Graphene\hspace{1mm}C_{6}$",'Interpreter','latex');
gp.set_xlabel("$X-Coordinates$",'Interpreter','latex');
gp.set_ylabel("$Y-Coordinates$",'Interpreter','latex');

%dumy = gp.draw("rect",5,5,3,4,'FaceColor','blue','Visible','off');
%a = gp.copy_to(dumy,10,10,10,2,'Visible','on');

cr = gp.draw("circle black",10,10,3);
copy_cr =  gp.copy_to(cr,20,10,'EdgeColor','black','FaceColor','red');

%rc = gp.draw("rect blue",13.5,18.5,3,3);
l = gp.draw("line",15,0,15,30,'Color','red');
bond = gp.draw("line black",0,0,0,0,'Visible','off');
lc = gp.copy_to(bond,0,20,30,20,'Visible','on');
gp.draw('crect blue',15,20,3,3);

%gp.set_text(l,'FBGENCER');

gp2 = lattice_drawer(figure(2),150,150);


%z = gp2.draw("cuboid",0,0,0,3,4,5);
%z2 = gp.draw("cuboid red",5,2,3,10,20,30,'FaceColor','blue');
for i = 10:10:100
	for j = 10
		p = gp2.draw("point yellow",i,j);
		t=gp2.set_text(p,sprintf("[%d,%d]",i,j));
	end
end
%v = gp2.draw("vector",0,0,30,30);
%ll = gp2.draw("line",10,10,10,20,20,20);

