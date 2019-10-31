classdef tightbinding < handle
   properties
      name = 0; %model name
      dim = 0; % dimension
      pvec = 0; %primitive vectors
      recip_vec = 0; %reciprocal vectors
      no_primvec = 0; %number of primitive vectors
      unit_cell = {};
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
       function tb = tightbinding(name,a1,varargin)
         tb.name = name;
         tb.set_primitive_vectors(a1,varargin{:});
       end
       %%
       function set_unit_cell(self,varargin)
       %nargin must be even, we will create a structure that holds our unit
       %cell-
            atom = struct('name','','pos',0);
            for i = 1:2:nargin-1
                atom.name = varargin{i};
                atom.pos = varargin{i+1};
               self.unit_cell{end+1} = atom;
            end
       end
       %%
        function set_primitive_vectors(self,a1,varargin)
          primvec = a1;
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
          disp(self.no_primvec);
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
        if(self.no_primvec > 1) 
          r.b2 = 2*pi*cross(a1,a3)/(dot(a2,cross(a1,a3)));
            if(self.no_primvec > 2)
              r.b3 = 2*pi*cross(a1,a2)/(dot(a3,cross(a1,a2)));
            end
        end
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
            
            bond.phase =  trans_vec;
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
                eig_matrix = zeros(size(self.unit_cell,2));
                if(size(eig_matrix) == 0)
                    error(message('Unit cell is undefined !'));
                end
                k = [kx(i,j),ky(i,j),kz(i,j)];
                for ob_iter = 1:bonding_no
                    bond = self.bonds{ob_iter};
                    q = exp(1i*dot(k,bond.phase*self.pvec));
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
            %plot(kx,self.E);
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
        
        precision = size(self.kvec{1},2);
        kx = [];
        ky = [];
        kz = [];
        no_point = nargin-3;
        for i = 1:no_point
            kx = [kx,linspace(varargin{i}(1),varargin{i+1}(1),precision)];
            ky = [ky,linspace(varargin{i}(2),varargin{i+1}(2),precision)];
            kz = [kz,linspace(varargin{i}(3),varargin{i+1}(3),precision)];
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
      function plot_lattice(self,gp,varargin)
            
            canvas_x = gp.canvas_x;
            canvas_y = gp.canvas_y;
            %gp.set_title('$Graphene\hspace{1mm}C_{6}$','Interpreter','latex');
            %gp.set_xlabel('$X$','Interpreter','latex');
            
            
            no_varg = nargin-2;
            if(no_varg>0)
                type_struct = varargin{1};
                plot_bond_obj = type_struct.bond;
                plot_atoms_obj = type_struct.atom;
            end
            
            cx = 0;
            cy = 0;
            gp.fig.CurrentAxes.XLim = [-canvas_x/2 canvas_x/2];
            gp.fig.CurrentAxes.YLim = [-canvas_y/2 canvas_y/2];
            
            normalize = 1e10;
            a1 = self.pvec(1,:).*normalize;
            if(self.no_primvec == 2) 
              a2 = self.pvec(2,:).*normalize;
            else
              a2 = zeros(1,size(a1,2));
            end

            
            %a3 = self.pvec(3,1);
            
            %atom_circles = {};
            %atom_bonds = {};
            
            %iterate for each atom in the unit cell and get the center
            %values
            for i = 1:size(self.unit_cell,2)
                atom = self.unit_cell{i};
                cx = cx+atom.pos(1)*normalize;
                cy = cy+atom.pos(2)*normalize;
            end
            
            main_cx = cx/size(self.unit_cell,2);
            main_cy = cy/size(self.unit_cell,2);
            
            lattice_fill_factorx = -canvas_x/2 : canvas_x/2;
            lattice_fill_factory = -canvas_y/2 : canvas_y/2 ;           
            
            %bunu 1d icin yaptim düzelmesi gerekiyor.
            if(a2 == 0)
                lattice_fill_factory = 1;
            end
            
            for e1 = lattice_fill_factorx
            for e2 = lattice_fill_factory
                cx = main_cx + e1 * a1(1) + e2 * a2(1);
                cy = main_cy + e1 * a1(2) + e2 * a2(2);
            for i = 1:size(self.bonds,2)
                bond = self.bonds{i};
                atom1 = self.unit_cell{bond.i};
                atom2 = self.unit_cell{bond.j};
                tvector = a1 .* bond.phase(1)+a2 .* bond.phase(2); % tvec translates our unitcell to new unitcells
                new_unitcell_cx = cx + tvector(1);
                new_unitcell_cy = cy + tvector(2);
                %gp.draw('vector',cx,cy,new_unitcell_cx,new_unitcell_cy);
                x1 = cx + atom1.pos(1) * normalize;
                y1 = cy + atom1.pos(2) * normalize;
                x2 = new_unitcell_cx + atom2.pos(1) * normalize;
                y2 = new_unitcell_cy + atom2.pos(2) * normalize;
                
                if(x1 == x2 && y1 == y2)
                   continue;
                end

                tbond = gp.copy_to(plot_bond_obj{1},x1,y1,x2,y2,'Visible','on');
                tsite1 = gp.copy_to(plot_atoms_obj{1},x1,y1,'Visible','on');
                tsite2 = gp.copy_to(plot_atoms_obj{2},x2,y2,'Visible','on');
            end
            end
            end
%             new_cx = cx;
%             new_cy = cy;
%             for i = 1
%             for j = 1
%                 new_cx = i*a1(1)+j*a2(1)+new_cx;
%                 new_cy = i*a1(2)+j*a2(2)+new_cy;
%                 for k = 1:size(atom_circles,2)
%                 %gp.draw('line black',cx-1,cy,cx+1,cy,'black');
%                     gp.copy_to(atom_circles{k},new_cx,new_cy);
%                 end
%                 new_cx = cx;
%                 new_cy = cy;
%             end
%             end

      end
   end
end


