% --------------------------------------------------------------------------------
%  MATLAB code used to generate the figures in the book:
%
%    "Quantum Transport: Atom to Transistor," by Supriyo Datta
%    published by Cambridge University Press, May 2005
%      (ISBN-10: 0521631459 | ISBN-13: 9780521631457)
%    http://www.cambridge.org/uk/catalogue/catalogue.asp?isbn=0521631459
%
% THIS FILE FOR: Chapter 5, Figure 5.4.1(a)
%
% --------------------------------------------------------------------------------
% Copyright (c) 2005  Supriyo Datta
% --------------------------------------------------------------------------------

clear all
%close all;

soa=.3787/3;soc=.0129/3;Esa=-8.3431;Epa=1.0414;Esc=-2.6569;Epc=3.6686;Esea=8.5914;
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

Ho=[h+h'	zeros(10);
     zeros(10)	h+h'];

hso=zeros(20);
hso(3,4)=-i*soa;
hso(3,15)=soa;
hso(4,15)=-i*soa;
hso(5,13)=-soa;
hso(5,14)=i*soa;
hso(6,7)=-i*soc;
hso(6,18)=soc;
hso(7,18)=-i*soc;
hso(8,16)=-soc;
hso(8,17)=i*soc;
hso(13,14)=i*soa;
hso(16,17)=i*soc;
Hso=hso+hso';

H = Ho+Hso;

eiglst = eig(H);
        E(Nk,:) = sort(real(eiglst));
		X(Nk)=-(Nk-1)/(Nt-1);%L-direction
		X1(Nk)=(Nk-1)/(Nt-1);%X-direction
end

hold on
if iter == 1
pl=plot(X,E);
else
pl=plot(X1,E);

end
axis([-1 1 -3 3])
xlabel('k (as fraction of maximum value)--->')
ylabel('Energy (eV) ---> ')
grid on
end