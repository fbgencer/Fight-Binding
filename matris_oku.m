f = fopen('matris_data.dat','r');


H = zeros(30);

for col = 1:30
	for row = 1:30
		s = fgetl(f)
		s = str2double(strsplit(s))
		H(row,col) =s(1)+1i*s(2); 
	end
end

fclose(f);