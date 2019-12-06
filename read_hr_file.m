function [bonds,deg_points] = read_hr_file(file_name)

f = fopen(file_name,'r');

%skip the comment line
fgetl(f);
matrix_row = str2double(fgetl(f)); % Number of wannier functions
nrpts = str2double(fgetl(f));

%read denegeracy points, number of degeneracy pts = nrpts

i = 1;
deg_points = zeros(1,nrpts);

while(i <= nrpts)
	s = fgetl(f);
	s = str2double(strsplit(s,' '));
	for j = 1:numel(s)
		if(isnan(s(j)) == 0)
			deg_points(i) = s(j);
			i = i + 1;
		end
	end 
end

%We know that matrix is square, we have matrix_row*matrix_row matrix and for each of them we will have an entry from
%hr file. So the size of bonds are already determined. But for each occurence of the same row and column we will make
%a matrix, for our tight binding class specification

% bonds = cell(matrix_row);

% for i = 1:numel(bonds)
% 	bonds{i} = {zeros(nrpts,3),zeros(nrpts,1)};
% end 

% bond = struct('phase',zeros(nrpts,3) ,'i',0,'j',0,'amp',zeros(nrpts,1));
% for i=1:numel(bonds)
% 	bonds{i} = bond;
% end

iter = 1;

tic;
bonds = textscan(f,'%f %f %f %d %d %f %f');
bonds{6} = bonds{6} + 1i*bonds{7};
bonds{7} = [];
%R1 = c{1};
%R2 = c{2};
%R3 = c{3};
%ar = c{6};
%ai = c{7};

%w = 1;

%Matlab holds in row order as the file hr.dat so instead holding m,n values we are putting them inside a cell at m,n
%1st one is R vector, second is amplitude.
% for i = 1:numel(ar)
% 	%yeni kod
% 	if(iter == numel(bonds)+1)
% 		iter = 1;
% 		w = w + 1;
% 	end
% 	bonds{iter}{1}(w,:) = [R1(i), R2(i), R3(i)];
% 	bonds{iter}{2}(w,1) = ar(i) + 1j*ai(i);	

% 	iter = iter + 1;

% 	% if(iter == numel(bonds)+1)
% 	% 	iter = 1;
% 	% 	w = w + 1;
% 	% end
% 	% bond = bonds{iter};

% 	% bond.phase(w,:) = [R1(i), R2(i), R3(i)];
% 	% bond.i = m(i);
% 	% bond.j = n(i);
% 	% bond.amp(w,1) = ar(i) + 1j*ai(i);	

% 	% bonds{iter} = bond;

% 	% iter = iter + 1;






% end
toc;

fclose(f);

end