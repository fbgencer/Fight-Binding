classdef tightbinding < handle
   properties
      name = 0; %model name
      dim = 0; % dimension
      pvec = 0; %primitive vectors
      recip_vec = 0; %reciprocal vectors
      no_primvec = 0; %number of primitive vectors
      hm = {}; % Hamiltonian matrix
      phase = {};
      no_unit = 0;
      kvec = {};
      E = 0; % without calculation it is just zero
      En = 0;
   end
   methods
       %selfect constructor
       %name: Name of the model
       %dim: Dimension of the model, 1d-2d-3d
       %unit_cell atom no
       %primvec: primitive vectors
       function tb = tightbinding(name,no_unit)
         tb.name = name;
         tb.no_unit = no_unit; % number of atoms in unit cell
       end
       %%
       function set_primitive_vectors(self,a1,varargin)
         if(nargin > 2)%include self,a1 and a2
            a2 = varargin{1};
            primvec = [a1;a2];
            if(nargin > 3 ) % include self,a1,a2,a3
                a3 = varargin{2};
                primvec = [a1;a2;a3];
            end
         end
         self.pvec = primvec;
         self.dim = size(primvec,2); % dimension of the lattice
         self.no_primvec = size(primvec,1); % number of primitive vectors 
       end
       %%
       function r = get_primitive_vectors(self)
        r.a1 = self.pvec(1,:);
        r.a2 = 0;
        r.a3 = 0;
        if(self.no_primvec > 1)
            r.a2 = self.pvec(2,:);
            if(self.no_primvec > 2)
                r.a3 = self.pvec(3,:);
            end
        end
       end
       %%
       function r = get_reciprocal_vectors(self)
        a1 = self.pvec(1,:);
        a2 = [0 1 0];
        a3 = [0 0 1];
        if(self.no_primvec > 1)
            a2 = self.pvec(2,:);
            if(self.no_primvec > 2)
                a3 = self.pvec(3,:);
            end
        end
        r.b1 = 2*pi*cross(a2,a3)/(dot(a1,cross(a2,a3)));
        r.b2 = 2*pi*cross(a1,a3)/(dot(a2,cross(a1,a3)));
        r.b3 = 2*pi*cross(a1,a2)/(dot(a3,cross(a1,a2)));
      end
       
%add_hopping(amplitude,interaction between pairs,translation vec)
%add_hopping(t,1,2,[0 0]);  % inside 0 A and B interaction interaction
%between zero and one which A and B respectively.
%add_hopping(t,1,2,[-1 0]);  % 0 -> 1
%add_hopping(t,2,1,[1 0]);   % 0 -> 2
%add_hopping(t,2,1,[0 1]);  % 0 -> 3
%add_hopping(t,1,2,[0 -1]);  % 0 -> 4
      
      function add_hopping(self,amp,index1,index2,trans_vec)
        if(index1 > 0 && index2 > 0)
            where = trans_vec * self.pvec;
            m = zeros(self.no_unit);
            m(index1,index2) =  amp;
            self.hm{end+1} = m;
            self.phase{end+1} = where;
        else
            disp('Atom indexes must be bigger than zero.');
        end
        
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function calculate_band(self)
        len = size(self.kvec{1},2);
        self.E = zeros(len,len);
        kx = self.kvec{1};
        ky = self.kvec{2};
        kz = self.kvec{3};
        bonding_no = size(self.hm,2);
        
        for i = 1:len
            for j = 1:len
                eig_matrix = zeros(self.no_unit);
                k = [kx(i,j),ky(i,j),kz(i,j)];
                for ob_iter = 1:bonding_no
                    m = self.hm{ob_iter};
                    v = self.phase{ob_iter};
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
                self.E(i,j) = eigens(1);
                %self.En(i,j) = eigens(2);
            end
        end
        %disp(eig_matrix);
        %disp(E);
        
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function set_kvector(self,from,to,len)
        k = linspace(from,to,len);
        if(self.no_primvec == 3)
            [kx,ky,kz] = meshgrid(k);
            self.kvec = {kx,ky,kz};
        elseif(self.no_primvec == 2)
            [kx,ky] = meshgrid(k);
            z = zeros(size(kx,1),size(kx,2));
            self.kvec = {kx,ky,z}; 
        elseif(self.no_primvec == 1)
            [kx] = meshgrid(k);
            z = zeros(size(kx,1),size(kx,2));
            self.kvec = {kx,z,z}; 
        end    
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function f = plot_band(self)
        kx = self.kvec{1};
        ky = self.kvec{2};
        kz = self.kvec{3};
        f = figure(1);
        if(self.no_primvec == 3)
            surf(kx,ky,kz,self.E);
        elseif(self.no_primvec == 2)
            surf(kx,ky,self.E);
            hold on;
            %surf(kx,ky,self.En);
        elseif(self.no_primvec == 1)
            surf(kx,self.E);
        end
        
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function f = plot_high_symmetry_points(self)
        a0 = 1.42*1e-10;
        a = a0*sqrt(3);
        f = figure(2);
        K1 = [-4*pi / (3*sqrt(3)*a), 0,0];
        Gamma = [0, 0, 0];
        M = [0, 2*pi / (3*a),0];
        K2 = [2*pi / (3*sqrt(3)*a), 2*pi / (3*a),0];
        
        len = size(self.kvec{1},2);
        kx = linspace(K1(1),Gamma(1),len);
        ky = zeros(1,len);
        kz = zeros(1,len);
       
        Energy = zeros(len,len);
        bonding_no = size(self.hm,2);
        
        for i = 1:len
            for j = 1:len
                eig_matrix = zeros(self.no_unit);
                k = [kx(i,j),ky(i,j),kz(i,j)];
                for ob_iter = 1:bonding_no
                    m = self.hm{ob_iter};
                    v = self.phase{ob_iter};
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
                Energy(i,j) = eigens(1);
                %self.En(i,j) = eigens(2);
            end
        end
        plot(kx,Energy);
      end
   end
end


