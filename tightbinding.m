classdef tightbinding < handle
   properties
      name = 0;
      dim = 0;
      pvec = 0;
      no_primvec = 0;
      hm = {}; % Hamiltonian matrix
      phase = {};
      no_unit = 0;
      kvec = {};
      E = 0; % without calculation it is just zero
   end
   methods
       %Object constructor
       %name: Name of the model
       %dim: Dimension of the model, 1d-2d-3d
       %unit_cell atom no
       %primvec: primitive vectors
       function tb = tightbinding(name,no_unit,primvec)
         tb.name = name;
         tb.pvec = primvec;
         tb.no_unit = no_unit; % number of atoms in unit cell
         tb.dim = size(primvec,2);
         tb.no_primvec = size(primvec,1); % number of primitive vectors
      end
%add_hopping(amplitude,interaction between pairs,translation vec)
%add_hopping(t,1,2,[0 0]);  % inside 0 A and B interaction interaction
%between zero and one which A and B respectively.
%add_hopping(t,1,2,[-1 0]);  % 0 -> 1
%add_hopping(t,2,1,[1 0]);   % 0 -> 2
%add_hopping(t,2,1,[0 1]);  % 0 -> 3
%add_hopping(t,1,2,[0 -1]);  % 0 -> 4
      
      function add_hopping(obj,amp,index1,index2,trans_vec)
        if(index1 > 0 && index2 > 0)
            where = trans_vec * obj.pvec;
            m = zeros(obj.no_unit);
            m(index1,index2) =  amp;
            obj.hm{end+1} = m;
            obj.phase{end+1} = where;
        else
            disp('Atom indexes must be bigger than zero.');
        end
        
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function calculate_band(obj)
        len = size(obj.kvec{1},2);
        obj.E = zeros(len,len);
        kx = obj.kvec{1};
        ky = obj.kvec{2};
        kz = obj.kvec{3};
        bonding_no = size(obj.hm,2);
        
        for i = 1:len
            for j = 1:len
                eig_matrix = zeros(obj.no_unit);
                k = [kx(i,j),ky(i,j),kz(i,j)];
                for ob_iter = 1:bonding_no
                    m = obj.hm{ob_iter};
                    v = obj.phase{ob_iter};
                    q = exp(1i*dot(k,v));
                    eig_matrix = eig_matrix + (m * q);
                    %fprintf("k = [%g %g %g] [%d]\n",k,ob_iter);
                    %disp(q);
                    %disp(v);
                    
                    %disp(eig_matrix);
                    %alpha = -t*(exp(-1i*dot(k,a1))+exp(-1i*dot(k,a2)) + 1);
                    %M = [Eo alpha; conj(alpha) Eo];
                end
                eigens = eig(eig_matrix);
                %disp(eigens);
                obj.E(i,j) = eigens(1);
            end
        end
        %disp(eig_matrix);
        %disp(E);
        
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function set_kvector(obj,from,to,len)
        k = linspace(from,to,len);
        if(obj.no_primvec == 3)
            [kx,ky,kz] = meshgrid(k);
            obj.kvec = {kx,ky,kz};
        elseif(obj.no_primvec == 2)
            [kx,ky] = meshgrid(k);
            z = zeros(size(kx,1),size(kx,2));
            obj.kvec = {kx,ky,z}; 
        elseif(obj.no_primvec == 1)
            [kx] = meshgrid(k);
            z = zeros(size(kx,1),size(kx,2));
            obj.kvec = {kx,z,z}; 
        end    
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function f = plot_band(obj)
        kx = obj.kvec{1};
        ky = obj.kvec{2};
        kz = obj.kvec{3};
        f = figure(1);
        if(obj.no_primvec == 3)
            surf(kx,ky,kz,obj.E);
        elseif(obj.no_primvec == 2)
            surf(kx,ky,obj.E);
        elseif(obj.no_primvec == 1)
            surf(kx,obj.E);
        end
        
        end
   end
end


