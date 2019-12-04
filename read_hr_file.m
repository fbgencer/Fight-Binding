function bonds = read_hr_file(file_name)

f = fopen(file_name,'r');

%skip the comment line
fgetl(f);
no_wannier = str2double(fgetl(f));
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


bonds = {};
while( ~ feof(f) )
	s = fscanf(f,'    %d    %d    %d   %d   %d    %f    %f',7);
	if(numel(s) ~= 7), continue, end

	bond.phase = [s(1) s(2) s(3)];
	bond.i = s(4);
	bond.j = s(5);
	bond.amp = s(6) + 1j*s(7);
	
	bonds{end+1} = bond;
end


%fscanf(f,'    %d    %d    %d    %d    %d    %f    %f',1)

fclose(f);

end