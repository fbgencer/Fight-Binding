clear; clc; close all;
%graphene consts
Eo = 1;
t = 2;
a0 = 1.42*1e-10;
a = a0*sqrt(3);

%primitive vectors
a1 = [a*sqrt(3)/2, a/2];
a2 = [a*sqrt(3)/2, -a/2];

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
%add_hopping(t,0,1,[0 0]);  % inside 0 A and B interaction interaction
%between zero and one which A and B respectively.
%add_hopping(t,0,1,[-1 0]);  % 0 -> 1
%add_hopping(t,1,0,[1 0]);   % 0 -> 2
%add_hopping(t,1,0,[0 1]);  % 0 -> 3
%add_hopping(t,0,1,[0 -1]);  % 0 -> 4


 


%set range in K-space
range = 2*pi/a;
len = 3;
k__ = linspace(-range,range,len);
[kx,ky] = meshgrid(k__);

E = zeros(len,len);

for i = 1:len
    for j = 1:len
        k = [kx(i,j),ky(i,j)];
        alpha = -t*(exp(-1i*dot(k,a1))+exp(-1i*dot(k,a2)) + 1);
        M = [Eo alpha; conj(alpha) Eo];
        eigens = eig(M);
        E(i,j) = eigens(2);
    end
end

f = figure(1);
surf(kx,ky,E);


figure(2);
plot(kx(1,:),E(1,:));
K1 = -4*pi / (3*sqrt(3)*a);
M = 2*pi / (3*a);
K2 = 2*pi / (3*sqrt(3)*a);
Gamma = 0;
xticks([K1 Gamma K2 M]);
xticklabels({'K','K','M','G'});
