function [bonds] = read_hr_file(file_name)

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

%hr File Format
%R(1) R(2) R(3) i j Amp_real Amp_imag


matrix_row2 = matrix_row*matrix_row;

R = zeros(nrpts,3);
matrix = zeros(matrix_row,matrix_row,nrpts);

tic;

c = textscan(f,'%f %f %f %d %d %f %f');
R1 = c{1};
R2 = c{2};
R3 = c{3};
m  = c{4};
n  = c{5};
ar = c{6};
ai = c{7};

iter = 0;

for i = 1:numel(c{1})

	if(mod(i,matrix_row2) == 1)
		%We have new R values completed NxN matrix
		iter = iter+1;
		R(iter,1) = R1(i);
		R(iter,2) = R2(i);
		R(iter,3) = R3(i);
	end

	matrix(m(i),n(i),iter) = complex(ar(i),ai(i))/deg_points(iter);

end

bonds.R = R;
bonds.matrix = matrix;


%bonds{6} = bonds{6} + 1i*bonds{7};
%bonds{7} = [];
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