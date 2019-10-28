classdef tightbinding < handle
   properties
      name = 0; %model name
      dim = 0; % dimension
      pvec = 0; %primitive vectors
      recip_vec = 0; %reciprocal vectors
      no_primvec = 0; %number of primitive vectors
      no_unit = 0;
      hm = {}; % Hamiltonian matrix
      phase = {};
      bonds = {};
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
       %%
%add_hopping(amplitude,interaction between pairs,translation vec)
%add_hopping(t,1,2,[0 0]);  % inside 0 A and B interaction interaction
%between zero and one which A and B respectively.
%add_hopping(t,1,2,[-1 0]);  % 0 -> 1
%add_hopping(t,2,1,[1 0]);   % 0 -> 2
%add_hopping(t,2,1,[0 1]);  % 0 -> 3
%add_hopping(t,1,2,[0 -1]);  % 0 -> 4
      
      function add_hopping(self,amp,index1,index2,trans_vec)
        if(index1 > 0 && index2 > 0)
            %Here we will hold a struct that contains i,j locations,
            %amplitude of the given hopping and phase value
            
            %where =
            %m = zeros(self.no_unit);
            %m(index1,index2) =  amp;
            %self.hm{end+1} = m;
            %self.phase{end+1} = where;
            bond.phase =  trans_vec * self.pvec;
            bond.i = index1;
            bond.j = index2;
            bond.amp = amp;
            self.bonds{end+1} = bond;
        else
            disp('Atom indexes must be bigger than zero.');
        end
        
      end
      %%
      function calculate_band(self)
          
        len = size(self.kvec{1},2);
        kx = self.kvec{1};
        ky = self.kvec{2};
        kz = self.kvec{3};
        
        [self.E,self.En] = calc_band_internal(self,kx,ky,kz);
        
      end
      function [Energy,negEnergy] = calc_band_internal(self,kx,ky,kz)
         %assume they have equal sizes
         ilen = size(kx,1);
         jlen = size(kx,2);
         Energy = zeros(ilen,jlen);
         negEnergy = Energy;
         bonding_no = size(self.bonds,2);
         for i = 1:ilen
            for j = 1:jlen
                eig_matrix = zeros(self.no_unit);
                k = [kx(i,j),ky(i,j),kz(i,j)];
                for ob_iter = 1:bonding_no
                    bond = self.bonds{ob_iter};
                    q = exp(1i*dot(k,bond.phase));
                    eig_matrix(bond.i,bond.j) = eig_matrix(bond.i,bond.j) +  (bond.amp * q);
                end
                eigens = eig(eig_matrix);
                Energy(i,j) = eigens(1);
                if(size(eigens,1) == 2)
                    negEnergy(i,j) = eigens(2);
                end
            end
         end
        
      end
      %%
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
      %%
      function f = plot_band(self,fig)
        kx = self.kvec{1};
        ky = self.kvec{2};
        kz = self.kvec{3};
        f = fig();
        s.positive_surface = 0;
        s.negative_surface = 0;
        if(self.no_primvec == 3)
            s.positive_surface = surf(kx,ky,kz,self.E);
        elseif(self.no_primvec == 2)
            s.positive_surface = surf(kx,ky,self.E);
            if(self.En ~= 0)
                hold on;
                s.negative_surface = surf(kx,ky,self.En);
            end
        elseif(self.no_primvec == 1)
            s.positive_surface = surf(kx,self.E);
        end
        f.UserData = s;
      end
      %%
      function f = plot_high_symmetry_points(self,fig,varargin)

        f = fig();
        
        if(nargin < 2)
            disp('In order to plot at least two high symmetry points are required.');
            return;
        end
        
        prec = size(self.kvec{1},2);
        kx = [];
        ky = [];
        kz = [];
        no_point = nargin-3;
        for i = 1:no_point
            kx = [kx,linspace(varargin{i}(1),varargin{i+1}(1),prec)];
            ky = [ky,linspace(varargin{i}(2),varargin{i+1}(2),prec)];
            kz = [kz,linspace(varargin{i}(3),varargin{i+1}(3),prec)];
        end
        
        [Energy,Energyn] = calc_band_internal(self,kx,ky,kz);
        p.positive_plot = plot(Energy);
        if(self.En ~= 0)
            hold on;
            p.negative_plot = plot(Energyn);
        end
        f.UserData = p;
      end
      %%
   end
end


