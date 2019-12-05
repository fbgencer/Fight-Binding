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

bonds = cell(matrix_row);
bond = struct('phase',[],'i',0,'j',0,'amp',[]);


for i=1:numel(bonds)
	bonds{i} = bond;
end

iter = 1;

w = 0;
while( ~ feof(f) )
	s = fscanf(f,'    %d    %d    %d   %d   %d    %f    %f',7);
	if(numel(s) ~= 7), continue, end
	
	if(iter == numel(bonds)+1)
		iter = 1;
	end

	bonds{iter}.phase(end+1,:) = [s(1) s(2) s(3)];
	bonds{iter}.i = s(4);
	bonds{iter}.j = s(5);
	bonds{iter}.amp(end+1,1) = s(6) + 1j*s(7);
	
	iter = iter + 1;

end


%fscanf(f,'    %d    %d    %d    %d    %d    %f    %f',1)

fclose(f);

end