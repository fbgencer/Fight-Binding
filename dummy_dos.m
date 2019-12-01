close all;
clear; clc;

step = 1;
kx = 0:step:100;
ky = kx;
kz = ky;
dim = 3;

if(dim == 1)
	ky = 0;
	kz = 0;
elseif(dim == 2)
	kz = 0;
end


E = zeros(1,size(kx,2)^dim);

ctr = 1;
for i1 = kx
	for i2 = ky
		for i3 = kz
			E(ctr)  = i1.^2+i2.^2+i3.^2;
			ctr = ctr + 1;
		end
	end
end
%E = kx.^2+ky.^2;

%plot(k,E)

h = histogram(E,500,'Visible','off');
counts = h.BinCounts;
x = linspace(h.BinLimits(1),h.BinLimits(2),size(counts,2)); 
i = 1;
% while( i <= numel(counts))
% 	if(counts(i) <= 1e-3)
% 		%disp("Evet ?")
% 		%disp(counts(i))
% 		counts(i) = [];
% 		x(i) = [];
% 		continue;
% 	end
% 	i = i+1;
% end
A = smoothdata(counts,'Gaussian');
f = figure(1);
A = smoothdata(A,'Gaussian');
plot(f.CurrentAxes,x,A,'-')
f.CurrentAxes.YLim(1) = 0;


% figure;

% E = sort(E);
% % E_sort = sort(E);


% allowed = [];
% occurance= [];

% for i = E;
% 	if(any(i == allowed))
% 		continue;
% 	end

% 	allowed(end+1) = i;
% 	occurance(end+1) = sum(E == i);
% end


% plot(allowed,occurance,'.-')