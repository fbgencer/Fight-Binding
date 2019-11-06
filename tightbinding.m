classdef tightbinding < handle
   properties
      name = 0; %model name
      dim = 0; % dimension
      pvec = 0; %primitive vectors
      recip_vec = 0; %reciprocal vectors
      no_primvec = 0; %number of primitive vectors
      unit_cell = {};
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
          primvec = [a1];
          if(nargin > 2)%include self,a1 and a2
            a2 = varargin{1};
            primvec = [a1;a2];
            if(nargin > 3 ) % include self,a1,a2,a3
                a3 = varargin{2};
                primvec = [a1;a2;a3];
            end
          end
          self.pvec = primvec;
          self.dim = size(primvec,1); % dimension of the lattice
          self.no_primvec = size(primvec,1); % number of primitive vectors 
          
        end
       %%
       function r = get_primitive_vectors(self)
        r{1} = self.pvec(1,:);
        r{2} = 0;
        r{3} = 0;
        if(self.no_primvec > 1)
            r{2} = self.pvec(2,:);
            if(self.no_primvec > 2)
                r{3} = self.pvec(3,:);
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
        r{1} = 2*pi*cross(a2,a3)/(dot(a1,cross(a2,a3)));
        if(self.no_primvec > 1) 
          r{2} = 2*pi*cross(a1,a3)/(dot(a2,cross(a1,a3)));
            if(self.no_primvec > 2)
              r{3} = 2*pi*cross(a1,a2)/(dot(a3,cross(a1,a2)));
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
      function calculate_band(self,kvec)
        
        kx = kvec{1};
        ky = kvec{2};
        kz = kvec{3};
        
        [self.E,self.En] = calc_band_internal(self,kx,ky,kz);
        
      end
      function [Energy,negEnergy] = calc_band_internal(self,kx,ky,kz)
        %assume they have equal sizes
        alen = size(kx,1);
        blen = size(kx,2);
        clen = size(kx,3);

        Energy = zeros(alen,blen,clen);
        negEnergy = Energy;
        bonding_no = size(self.bonds,2);
        if(size(self.unit_cell,2) == 0)
          error(message('Unit cell is undefined !'));
        end

         for a = 1:alen
            for b = 1:blen
              for c = 1:clen
                  eig_matrix = zeros(size(self.unit_cell,2));
                  k = [kx(a,b,c),ky(a,b,c),kz(a,b,c)];
                  for ob_iter = 1:bonding_no
                      bond = self.bonds{ob_iter};
                      q = exp(1i*dot(k,bond.phase*self.pvec));

                      eig_matrix(bond.i,bond.j) = eig_matrix(bond.i,bond.j) +  (bond.amp * q);
                  end
                  eigens = eig(eig_matrix);
                  Energy(a,b,c) = real(eigens(1));

                  if(size(eigens,1) == 2)
                    negEnergy(a,b,c) = real(eigens(2));
                  end
                end
            end
         end
        
      end
      %%
      function kvec = set_kvector(self,from,to,len)
        k = linspace(from,to,len);
        if(self.no_primvec == 3)
            [kx,ky,kz] = meshgrid(k);
            kvec = {kx,ky,kz};
        elseif(self.no_primvec == 2)
            [kx,ky] = meshgrid(k);
            z = zeros(size(kx,1),size(kx,2));
            kvec = {kx,ky,z}; 
        elseif(self.no_primvec == 1)
            [kx] = meshgrid(k);
            z = zeros(size(kx,1),size(kx,2));
            kvec = {kx,z,z}; 
        end    
      end
      %%
      function s = plot_energy_band(self,fig,kvec,plot_type,varargin)

        kx = kvec{1};
        ky = kvec{2};
        kz = kvec{3};

        [Ep,En] = calc_band_internal(self,kx,ky,kz);
        f = fig();
        s.positive_surface = 0;
        s.negative_surface = 0;

        if(plot_type == 'surface')  
          if(self.no_primvec == 3)
              s.positive_surface = surf(kx(:,:,1),ky(:,:,1),Ep(:,:,1),varargin{:});
              if(En ~= 0)
                  hold on;
                  s.negative_surface = surf(kx(:,:,1),ky(:,:,1),En(:,:,1),varargin{:});
              end
              
          elseif(self.no_primvec == 2)
              s.positive_surface = surf(kx,ky,Ep,varargin{:});
              if(En ~= 0)
                  hold on;
                  s.negative_surface = surf(kx,ky,En,varargin{:});
              end
          elseif(self.no_primvec == 1)
              s.positive_surface = surf(kx,Ep,varargin{:});
          else
            error(message('Before plotting, define primitive vectors'));
          end
          f.UserData = s;

        elseif(plot_type == 'contour')
          if(self.no_primvec == 3)
              s.positive_surface = surf(kx(:,:,1),ky(:,:,1),Ep(:,:,1),varargin{:});
              
          elseif(self.no_primvec == 2)
              s.positive_contour = contour(kx,ky,Ep,varargin{:});
              if(En ~= 0)
                hold on;
                s.positive_contour = contour(kx,ky,En,varargin{:});
              end
          elseif(self.no_primvec == 1)
              s.positive_contour = contour(kx,Ep,varargin{:});
          else
            error(message('Before plotting, define primitive vectors'));
          end
          f.UserData = s;          
        end

      end
      %%
      function f = plot_high_symmetry_points(self,fig,varargin)

        f = fig();
        
        if(nargin < 2)
            disp('In order to plot at least two high symmetry points are required.');
            return;
        end
        precision = 100;
        kx = [];
        ky = [];
        kz = [];
        no_point = nargin-3;
        for i = 1:no_point
            kx = [kx,linspace(varargin{i}(1),varargin{i+1}(1),precision)];
            ky = [ky,linspace(varargin{i}(2),varargin{i+1}(2),precision)];
            kz = [kz,linspace(varargin{i}(3),varargin{i+1}(3),precision)];
        end

        [Ep,En] = calc_band_internal(self,kx,ky,kz);
        p.positive_plot = plot(Ep);
        if(En ~= 0)
            hold on;
            p.negative_plot = plot(En);
        end
        f.UserData = p;
      end
      %%
      function plot_only_atoms(self,gp,varargin)
        x_offset = 0;
        y_offset = 0;
        z_offset = 0;

        normalize = 1e10;
        gp.xaxis_symmetric();
        gp.yaxis_symmetric();
        gp.zaxis_symmetric();
        

        if(nargin > 2)
          atom_object = varargin{1}.atoms;
        else
          atom_object = cell(1,size(self.unit_cell,2));
          colors = {"red","blue","green","magenta","yellow","black"};
          for i = 1:size(atom_object,2)
            atom_object{i} = gp.draw("circle "+colors{mod(i,6)+1},0,0,0.2,'Visible','off');
          end
        end



        %first get the atoms from the unitcell
        cx = 0; % unitcell center value
        cy = 0;
        cz = 0;
        for i = 1:size(self.unit_cell,2)
            atom = self.unit_cell{i}; %unit_cell{} contains atoms inside the unit cell
            cx = cx + atom.pos(1) .* normalize;
            cy = cy + atom.pos(2) .* normalize;
            if(size(atom.pos,2) > 2)
              cz = cz + atom.pos(3) .* normalize;
              
            end
        end
        

        for i = 1:size(self.unit_cell,2)
          atom = self.unit_cell{i};
          x = x_offset + cx + atom.pos(1)*normalize;
          y = y_offset + cy + atom.pos(2)*normalize;
          z = 0;

          if(size(atom.pos,2) > 2)
            z = z_offset + cz + atom.pos(3) * normalize;
          end

          if(size(atom.pos,2) > 2)
            gp.copy_to(atom_object{i},x,y,z,'Visible','on');
          else 
            gp.copy_to(atom_object{i},x,y,'Visible','on');
          end          
            
        end

      end
      %%
      function plot_only_bonds(self,gp,varargin)
        x_offset = 0;
        y_offset = 0;
        z_offset = 0;
        normalize = 1e10;
        gp.xaxis_symmetric();
        gp.yaxis_symmetric();
        gp.zaxis_symmetric();
        
        
        normalize = 1e10;
        a1 = self.pvec(1,:).*normalize;
        if(self.no_primvec > 1) 
          a2 = self.pvec(2,:).*normalize;
          if(self.no_primvec > 2)
            a3 = self.pvec(3,:).*normalize;
          else
            a3 = zeros(1,size(a1,2));
          end
        else
          a2 = zeros(1,size(a1,2));
        end
        if(nargin > 2)
          line_object = varargin{1}.bonds{1};
        else
          line_object = gp.draw("line black",0,0,0,0,'Visible','off');
        end

        for i = 1:size(self.bonds,2)
          bond = self.bonds{i};
          
          if(bond.phase == 0 & bond.i == bond.j)
            continue;
          end

          % tvec translates our unitcell to new unitcells
          tvector = a1 .* bond.phase(1);
          if(size(bond.phase,2) > 1)
            tvector = tvector + a2 .* bond.phase(2);
            if(size(bond.phase,2) > 2)
              tvector = tvector + a3 .* bond.phase(3);
            end
          end 

          atom1 = self.unit_cell{bond.i};
          atom2 = self.unit_cell{bond.j};
          x(1) = x_offset + atom1.pos(1) .* normalize;
          y(1) = y_offset + atom1.pos(2) .* normalize;
          x(2) = x_offset + tvector(1) + atom2.pos(1) .* normalize;
          y(2) = y_offset + tvector(2) + atom2.pos(2) .* normalize;
          z = [0 0];

          if(size(atom1.pos,2) > 2)
            z(1) = z_offset + atom1.pos(3) .* normalize;
            z(2) = z_offset + tvector(3) + atom2.pos(3) .* normalize;
          end

          if(x(1) == x(2) & y(1) == y(2) &  z(1) == z(2))
            continue;
          end

          if(x > gp.canvas_x/2 | x < -gp.canvas_x/2 | y > gp.canvas_y/2 | y < -gp.canvas_y/2)
            continue;
          end

          fprintf( "[%g %g][%g %g][%g %g]\n",x,y,z );

          gp.copy_to(line_object,x(1),y(1),z(1),x(2),y(2),z(2),'Visible','on');

        end

      end
      %%

      function plot_lattice(self,gp,varargin)

        x_offset = 0;
        y_offset = 0;
        z_offset = 0;

        line_object = "None";
        atom_object = "None";

        lattice_fill_factorx = -gp.canvas_x/2 : gp.canvas_x/2;
        lattice_fill_factory = -gp.canvas_y/2 : gp.canvas_y/2;
        lattice_fill_factorz = -gp.canvas_z/2 : gp.canvas_z/2;
        gp.xaxis_symmetric();
        gp.yaxis_symmetric();
        gp.zaxis_symmetric();

        i = 1;
        while(i <= size(varargin,2))
          st = string(varargin{i});
          i = i+1;
          switch st
            case 'atoms'
              atom_object = varargin{i};
            case 'bonds'
              line_object = varargin{i}{1};
            case 'x'
              lattice_fill_factorx = varargin{i};
              %gp.change_axis_limits("x",lattice_fill_factorx(1),lattice_fill_factorx(end));
            case 'y'
              lattice_fill_factory = varargin{i};
              %gp.change_axis_limits("y",lattice_fill_factory(1),lattice_fill_factory(end));
            case 'z'
              lattice_fill_factorz = varargin{i}; 
              %gp.change_axis_limits("z",lattice_fill_factorz(1),lattice_fill_factorz(end));        
            otherwise
              disp('eRRor');
              return;
          end
          i = i+1;
        end

        if(isstring(line_object) & line_object == "None");
          line_object = gp.draw("line black",0,0,0,0,'Visible','off');
        end
        if(isstring(atom_object) & atom_object == "None")
          atom_object = cell(1,size(self.unit_cell,2));
          colors = {"red","blue","green","magenta","yellow","black"};
          for i = 1:size(self.unit_cell,2)
            atom_object{i} = gp.draw("point "+colors{mod(i,6)+1},0,0,0.2,'Visible','off');
          end
        end

        if(self.dim == 1)
          lattice_fill_factory = 1;
        end


        
        normalize = 1e10;
        a1 = self.pvec(1,:).*normalize;
        if(self.no_primvec > 1) 
          a2 = self.pvec(2,:).*normalize;
        else
          a2 = zeros(1,size(a1,2));
        end
        if(self.no_primvec > 2)
          a3 = self.pvec(3,:).*normalize;
        else
          a3 = zeros(1,size(a1,2));
        end


        cx = 0;
        cy = 0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BOND plottign
        old_x = 0;
        old_y = 0;
        old_z = 0;
        for i = 1:size(self.bonds,2)
          for e1 = lattice_fill_factorx
          for e2 = lattice_fill_factory
          for e3 = lattice_fill_factorz
            tx = e1 * a1(1) + e2 * a2(1) + e3 * a3(1);
            ty = e1 * a1(2) + e2 * a2(2) + e3 * a3(2);
            tz = e1 * a1(3) + e2 * a2(3) + e3 * a3(3);
          
            bond = self.bonds{i};
            
            if(bond.phase == 0 & bond.i == bond.j)
              continue;
            end

            % tvec translates our unitcell to new unitcells
            tvector = a1 .* bond.phase(1);
            if(size(bond.phase,2) > 1)
              tvector = tvector + a2 .* bond.phase(2);
              if(size(bond.phase,2) > 2)
                tvector = tvector + a3 .* bond.phase(3);
              end
            end 
            atom1 = self.unit_cell{bond.i};
            atom2 = self.unit_cell{bond.j};
            x(1) = x_offset + atom1.pos(1) * normalize;
            y(1) = y_offset + atom1.pos(2) * normalize;
            x(2) = x_offset + tvector(1) + atom2.pos(1) * normalize;
            y(2) = y_offset + tvector(2) + atom2.pos(2) * normalize;
            x = x+tx;
            y = y+ty;
            z = [0 0];

            if(size(atom1.pos,2) > 2)
              z(1) = z_offset + atom1.pos(3) * normalize;
              z(2) = z_offset + tvector(3) + atom2.pos(3) * normalize;
              z = z+tz;
            end

            
            if(x == old_x & y == old_y & z == old_z)
              continue;
            end
            if(x(1) == x(2) & y(1) == y(2) &  z(1) == z(2))
              continue;
            end

            if(x < lattice_fill_factorx(1) | x > lattice_fill_factorx(end) | y < lattice_fill_factory(1) | y > lattice_fill_factory(end) | ...
               z > lattice_fill_factorz(1) | z < lattice_fill_factorz(end) )
              continue;
            end

            gp.copy_to(line_object,x(1),y(1),z(1),x(2),y(2),z(2),'Visible','on');

            old_x = x;
            old_y = y;
            old_z = z;
          end
          end
          end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %first get the atoms from the unitcell
        cx = 0; % unitcell center value
        cy = 0;
        for i = 1:size(self.unit_cell,2)
            atom = self.unit_cell{i}; %unit_cell{} contains atoms inside the unit cell
            cx = cx+atom.pos(1)*normalize;
            cy = cy+atom.pos(2)*normalize;
        end

        old_x = 0;
        old_y = 0;
        old_z = 0;

        for i = 1:size(self.unit_cell,2)
          for e1 = lattice_fill_factorx
            for e2 = lattice_fill_factory
              for e3 = lattice_fill_factorz

                tx = e1 * a1(1) + e2 * a2(1) + e3 * a3(1);
                ty = e1 * a1(2) + e2 * a2(2) + e3 * a3(2);
                tz = e1 * a1(3) + e2 * a2(3) + e3 * a3(3);

                atom = self.unit_cell{i};
                x = tx + x_offset + cx + atom.pos(1)*normalize;
                y = ty + y_offset + cy + atom.pos(2)*normalize;
                z = 0;

                if(size(atom.pos,2) > 2)
                  z = atom.pos(3) * normalize;
                  z = z+tz;
                end

                if(x == old_x & y == old_y & z == old_z)
                  continue;
                end

                if(x < lattice_fill_factorx(1) | x > lattice_fill_factorx(end) | y < lattice_fill_factory(1) | y > lattice_fill_factory(end) | ...
                   z > lattice_fill_factorz(1) | z < lattice_fill_factorz(end) )
                  continue;
                end

                if(size(atom.pos,2) > 2)
                  gp.copy_to(atom_object{i},x,y,z,'Visible','on');
                else 
                  gp.copy_to(atom_object{i},x,y,'Visible','on');
                end
                
                old_x = x;
                old_y = y;
                old_z = z;
              end
            end
          end
        end
      end
        %%
      function ln = plot_brillouin_zone(self,gp,varargin)
        %ln is line of the brillouin zone, so user can change its properties later

        plot_points = 0;
        plot_lines = 0;
        plot_coordinates = 0;

        nvargin = size(varargin,2);
        for i = 1:nvargin
          if(varargin{i} == "plot points")
            plot_points = 1;
            if(i < nvargin & ischar(varargin{i+1}) == 0 & isstring(varargin{i+1}) == 0)
              i = i+1;
              point_obj = varargin{i};
            else
              point_obj = gp.draw('point black',0,0,0.5,'Visible','off');
            end
          elseif(varargin{i} == "plot lines")
            plot_lines = 1;
            if(i < nvargin & ischar(varargin{i+1}) == 0 & isstring(varargin{i+1}) == 0)
              i = i+1;
              line_obj = varargin{i};
            else
              line_obj = gp.draw('line red',0,0,0,0,'Visible','off','ZData',0.5,'LineWidth',2);
            end
          elseif(varargin{i} == "plot coordinates")
            plot_coordinates = 1;          
          end
        end


        b = self.get_reciprocal_vectors();
        %Create closest points using recip vectors
        %assumption is we can create brillouin zone with the -1,0,1 indices of recp vectors
        points = [];
        if(self.dim == 2)
          for i = [-1,0,1]
            for j = [-1,0,1]
              p = i.*b{1}+j.*b{2};
              points = [points, [p(1);p(2)] ];
            end
          end
    
          %using voronoi calculate region vertices
          [vpx,vpy] = voronoi(points(1,:),points(2,:));
          %Voronoi calculates a lot of points but chose the first (XXXTODO think about this later, is it ok ?)
          vpx1 = vpx(1,:);
          vpy1 = vpy(1,:);

          %now we will eliminate same vertex points
          i = 1;
          while( i <= size(vpx1,2) )
            pt(1) = vpx1(i);
            pt(2) = vpy1(i);
            j = i+1;
            while(j <= size(vpx1,2) )
              if(pt(1) == vpx1(j) & pt(2) == vpy1(j))
                vpx1(j) = [];
                vpy1(j) = [];
              end
              j = j+1;
            end
            i = i+1;
          end
          %Now vpx1 and vpy1 holds different vertex points, which one is closes to the origin ?, using hypot we can find closes points

          hypots = transpose(hypot(vpx1(:),vpy1(:)));
          %now sort hypots and sort vpx1 and vpy1 according to hypot
          [hypots,I] = sort(hypots);
          vpx1 = vpx1(I);
          vpy1 = vpy1(I);

          %get hypot minimum
          hypotmin = min(hypots);
          i = 1;
          while (i <= size(hypots,2))
            if(abs(hypots(i)-hypotmin) > 1e-4) % k-space is huge and truncation will cause some equality errors,our tolerance is 1e-4, 
              %if there is a difference between hypot values then we will eliminate them until finding the closest points
              hypots(i) = [];
              vpx1(i) = [];
              vpy1(i) = [];
            else
              i = i + 1;
            end 
          end

          %now vpx1 and vpy1 hold closes points but in order to draw the brillouin zone we need to sort them in terms of angle between origin
          angles = atan2d(vpy1,vpx1);
          [angles,I] = sort(angles);
          vpx1 = vpx1(I);
          vpy1 = vpy1(I);

          if(plot_points)
            for i = 1:size(vpx1,2)
              gp.copy_to(point_obj,vpx1(i),vpy1(i),'Visible','on');
            end
          end

          if(plot_coordinates)
            for i = 1:size(vpx1,2)
              x = fix(vpx1(i));
              y = fix(vpy1(i));
              txt = pp(x,2) + "," + pp(y,2);
              %txt = sprintf("$$%g,%g$$",x,y);
              

              %txt = replace(txt,",","},");
              %txt = replace(txt,"]","}]");
              t = gp.draw('text',x,y,txt,'Color','black','Interpreter','Latex');
              t.Position(3) = 0.5;
            end            
          end

          if(plot_lines)
            %finally connect first value to the array
            vpx1(end+1) = vpx1(1);
            vpy1(end+1) = vpy1(1);
            %and now plot as a line
            ln = gp.copy_to(line_obj,vpx1,vpy1,'Visible','on');
          end

        end

      end
   end
end

function s=pp(x,n)
% pretty-print a value in scientific notation
% usage: s=pp(x,n)
% where: x a floating-point value
%        n is the number of decimal places desired
%        s is the string representation of x with n decimal places, and
%          exponent k
if(x ~= 0)
  exponent=floor(log10(abs(x))); %to accomodate for negative values
  mantissa=x/(10^exponent);
  s=sprintf('$%*.*f \\times 10^{%d}$',n+3,n,mantissa,exponent);
else
  s = sprintf('0');
end
% returns something like '$1.42 \times 10^{-1}$'
end