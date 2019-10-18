clear; clc; close all;
%graphene consts
Eo = 1;
t = 2;
a0 = 1.42*1e-10;
a = a0*sqrt(3);

%primitive vectors
a1 = [a*sqrt(3)/2, a/2,  0];
a2 = [a*sqrt(3)/2, -a/2, 0];

tb = tightbinding('graphene',2,[a1;a2]);


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

tb.add_hopping(Eo,1,1,[0 0]);  % 0 to 0
tb.add_hopping(Eo,2,2,[0 0]);  % 0 to 0

tb.add_hopping(-t,1,2,[0 0]);  % 0 to 0
tb.add_hopping(-t,2,1,[0 0]);  % 0 to 0
tb.add_hopping(-t,1,2,[-1 0]);  % 0 to 1 with -a1 vector
tb.add_hopping(-t,2,1,[1 0]);   % 0 to 2 with a1 vector
tb.add_hopping(-t,1,2,[0 -1]);  % 0 to 3
tb.add_hopping(-t,2,1,[0 1]);  % 0 to 4


range = 5*pi/a;
len = 80;
tb.set_kvector(-range,range,len);
tb.calculate_band();
f = tb.plot_band();
%f.Color = [0.5 0.5 0.2];
