clear all;
clc;

A1 = sym('A1');
A2 = sym('A2');
C1 = sym('C1');
C2 = sym('C2');
B = sym('B');

get_diag = @(A) A(sub2ind(size(A),1:size(A,1),1:size(A,2))); 

H = [A1,B,C1,0; conj(B), A2, 0, C2; C1, 0, -A1, -B; 0, C2, -conj(B), -A2];
[eigvec,eigval] = eig(H);
eigval = get_diag(eigval);


V1 = eigvec(:,1);
V2 = eigvec(:,2);
V3 = eigvec(:,3);
V4 = eigvec(:,4);

L1 = eigval(1);
L2 = eigval(2);
L3 = eigval(3);
L4 = eigval(4);


%HH = [A1 A2; C1 C2];
%[hvec,heig] = eig(HH);
%heig = get_diag(heig);