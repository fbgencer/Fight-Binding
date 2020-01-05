classdef tightbinding < handle
   properties
      dim = 0; % dimension
      orbitals = {};
      no_orbital = 1;
      pvec = 0; %primitive vectors
      recip_vec = 0; %reciprocal vectors
      no_primvec = 0; %number of primitive vectors
      unit_cell = {};
      phase = {};
      bonds = {};
      kvec = {};
      E = {}; % without calculation it is just zero

      %% Spin orbit coupling - SOC
      is_soc = 0; %we will open it after user calls add_soc

      %normalization unit for drawing lattices and k space
      spatial_unit = 1;
      spatial_unit_str = 'nm';

      deg_points = [];

      Efermi = 0;


      hr_file = '';

   end
   methods
       %selfect constructor
       %name: Name of the model
       %dim: Dimension of the model, 1d-2d-3d
       %unit_cell atom no
       %primvec: primitive vectors
       function tb = tightbinding(dim,a1,varargin)
          if(isnumeric(dim) == 0)
            error('Problem dimension must be numeric')
          end
          dim = int8(dim); %convert to integer..
          if(dim > 3)
            error('Dimension cannot be bigger than 3')
          end
          tb.dim = dim;
          tb.set_primitive_vectors(a1,varargin{:});
       end
       %%
       function set_unit_cell(self,varargin)
       %nargin must be even, we will create a structure that holds our unit
       %cell- 
            atom = struct('name','','pos',0);
            for i = 1:2:nargin-1
                atom.name = varargin{i};
                atom.pos = varargin{i+1}*self.spatial_unit;
               self.unit_cell{end+1} = atom;
            end
       end
       %%
       function set_orbital(self,varargin)

          self.no_orbital = 0;

          for i = 1:size(varargin,2)
            if(isstring(varargin{i}) | ischar(varargin{i}))
              orb_cell = strsplit(varargin{i},',');
              self.orbitals{end+1} = orb_cell;
              self.no_orbital = self.no_orbital + numel(orb_cell);
            else
              error('Orbital names must be string');
            end
          end
          %self.no_orbital = size(self.orbitals,2);
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
          self.no_primvec = size(primvec,1); % number of primitive vectors 
          
        end
       %%
       function set_fermi_level(self,ef)
          self.Efermi = ef;
       end
       %%
       function set_metric_unit(self,m)
        
        switch(m)
           case 'm'
              self.spatial_unit = 1; 
               self.spatial_unit_str = m;
           case 'nm'
               self.spatial_unit = 1e-9;
                self.spatial_unit_str = m;
           case 'A'
               self.spatial_unit = 1e-10;
                self.spatial_unit_str = m;
               
            otherwise
                error('Undefined unit, available units are {meter,nanometer,angstrom}');
        end
        self.pvec = self.pvec .* self.spatial_unit;
        for i = 1:size(self.unit_cell,2)
            self.unit_cell{i}.pos = self.unit_cell{i}.pos * self.spatial_unit;
        end
       
       
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
      function spin_orbit_coupling(self,state)
        if(strcmp('true',state) | state == 1)
          self.is_soc = 1;
        elseif(strcmp('false',state) | state == 0)
          self.is_soc = 0;
        end
      end
       %%
       function add_hrfile(self,filename,varargin)
        %Do not try to read again and again unless user specifically makes self.bonds = {}
        if(isempty(self.bonds))
          disp('Reading hr file...');
          tic;
            self.hr_file = filename;
            self.bonds = read_hr_file(filename);
            toc;
          disp('hr file is read succesfully.')
        end
       end
       %%
      function add_hopping(self,amp,index1,index2,trans_vec,varargin)
        %varargin part for orbitals, size must be 2




        spin1 = 0; %means undefined;
        spin2 = 0; %means undefined;

        % If user did not specify spin we will assume that term includes both spin terms, in the matrix we will write for both spin terms
        % But for specific spin states it is gonna write to exact row-column
        %check for spin up or down
        %This looks so bad, we need to change..XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        if(self.is_soc)
          where = strfind(index1,'+');
          if(where > 1) 
            index1 = index1(1:where-1);
            spin1 = 1;
          else
            where = strfind(index1,'-');
            if(where > 1)
              spin1 = -1;
              index1 = index1(1:where-1);
            end
          end
          where = strfind(index2,'+');
          if(where > 1) 
            spin2 = 1;
            index2 = index2(1:where-1);
          else
            where = strfind(index2,'-');
            if(where > 1)
              spin2 = -1;
              index2 = index2(1:where-1);
            end
          end              
        end


        %convert string indexes to numerical values
        if( isstring(index1) | ischar(index1) | isstring(index2) | ischar(index2))
          for i = 1:size(self.unit_cell,2)
            atom = self.unit_cell{i};
            if(strcmp(index1,atom.name))
              %bond.atoms{1} = i;
              index1 = i;
            end
            if(strcmp(index2,atom.name))
              %bond.atoms{2} = i;
              index2 = i;
            end
          end
        end

        if( isstring(index1) | ischar(index1) | isstring(index2) | ischar(index2)  )
          error('Undefined atom name');
        end   

        bond.atoms = {index1,index2};



        where = find(strcmp(varargin,'mode'));
        overwrite = 0; append = 0;
        if(where)
          overwrite = find(strcmp(varargin,'overwrite'));
          append = find(strcmp(varargin,'append'));
          if(overwrite == where+1), overwrite = 1;
          elseif(append == where+1), append = 1;
          else, error('Undefined mode. Append or overwrite is available'); end
        end

        where = find(strcmp(varargin,'sym'));
        if(where)
          symbolic_amp = varargin{where+1};
          if(isstring(symbolic_amp) | ischar(symbolic_amp))
            symbolic_amp = str2sym(symbolic_amp);
          end
          if(where == 1)
            where = 3;
          elseif(where == 3)
            where = 1;
          end
        else
          symbolic_amp = 'None';
          where = 1;
        end

  
        bond.orbitals = {'s','s'};

        offset1 = 0;
        offset2 = 0;
        orbital1 = 0;
        orbital2 = 0;

        %Check whether given is orbital string or symbolic entry
        if(size(varargin,2) > 2 | (size(varargin,2) == 2 & symbolic_amp == 0))
          
          for i = 1:size(self.orbitals,2)
            for orb_iter = 1:numel(self.orbitals{i})
              %disp(varargin{where});disp(self.orbitals{i}{orb_iter})
                if(orbital1 == 0)
                  if( strcmp(varargin{where},self.orbitals{i}{orb_iter}) )
                    bond.orbitals{1} = varargin{where};
                    orbital1 = i;
                    index1 = offset1 + numel(self.orbitals{i}) .* (index1-1) + orb_iter;
                  elseif(orb_iter == numel(self.orbitals{i}))
                      offset1 = offset1 + size(self.unit_cell,2) * numel(self.orbitals{i});
                  end
                end
                if(orbital2 == 0)
                  if( strcmp(varargin{where+1},self.orbitals{i}{orb_iter}) )
                    bond.orbitals{2} = varargin{where+1};
                    orbital2 = i;
                    index2 = offset2 + numel(self.orbitals{i}) .* (index2-1) + orb_iter;
                  elseif(orb_iter == numel(self.orbitals{i}))
                      offset2 = offset2 + size(self.unit_cell,2) * numel(self.orbitals{i});
                  end
                end
            end            
          end
        end
        
        %Index shifting due to soc
        %We make indexes as matrices
        if(self.is_soc)
          matrix_row_size = size(self.unit_cell,2)*self.no_orbital; %both spin up + spin down states
          if(spin1 == 0)
            index1 = [index1, index1 + matrix_row_size];
          end
          if(spin2 == 0)
            index2 = [index2, index2 + matrix_row_size];
          end
          if(spin1 < 0), index1 = index1 + matrix_row_size; end
          if(spin2 < 0), index2 = index2 + matrix_row_size; end
        end

        %index1
        %index2

        for i1 = 1:numel(index1)
          if(any(index1) > 0 && any(index2 > 0) )
              %Here we will hold a struct that contains i,j locations,
              %amplitude of the given hopping and phase value

              %First check the translation vector (trans_vec) and number of primitive vectors
              if(size(trans_vec,2) ~= size(self.pvec,1) )
                error('Problem dimension and translation vector size does not match!');
              end

              multiple_entry_flag = (size(amp,1) > 1 & size(amp,1) == size(trans_vec,1));
              
              ientry=1;
              %for ientry = 1:size(amp,1)
                where = numel(self.bonds)+1;

                %Also check that maybe user is trying add same elements, we should not add more bonds, we will quietly return
                for iter = 1:numel(self.bonds)
                  temp_bond = self.bonds{iter};
                  if(temp_bond.i == index1(i1) & temp_bond.j == index2(i1) & temp_bond.phase == trans_vec(ientry,:))
                    %disp(any(append));
                    %disp(iter)
                    if(append == 0 && overwrite == 0)
                      continue;
                    end
                    if(any(append) == 1)
                      amp = self.bonds{iter}.amp + amp;
                      %amp
                      if(symbolic_amp == 'None'), symbolic_amp = self.bonds{iter}.symbolic_amp;
                      else, symbolic_amp = self.bonds{iter}.symbolic_amp + symbolic_amp; end

                      where = iter;
                      %where
                    end
                  end
                end

                %Here we will create hermitian matrix so without expecting from user, we can add hermitian conjugate here
                bond.phase = trans_vec;%trans_vec(ientry,:);
                bond.i = index1(i1);
                bond.j = index2(i1);
                bond.amp = amp;%amp(ientry,1);
                bond.symbolic_amp = symbolic_amp;
                self.bonds{where} = bond;
                where = where + 1;

                %disp('Added');
                %bond

                %fprintf('index %d %d\n',index1,index2 );
                %We can give user some options to close this flag, LATER ADD this
                flag_create_hermitian_conj = 1;
                %Add only Non-Diagonal elements
                if((all(index1 ~= index2) | any(bond.phase ~= -trans_vec)) )
                  bond.phase = -trans_vec;%trans_vec(ientry,:);
                  bond.i = index2(i1);
                  bond.j = index1(i1);
                  bond.amp = conj(amp);%conj(amp(ientry,1));
                  if(isOctave == 0)
                    bond.symbolic_amp = conj(symbolic_amp);
                  end
                  %now change the atoms as well
                  y = bond.atoms;
                  bond.atoms{1} = y{2};
                  bond.atoms{2} = y{1}; 
                  self.bonds{where} = bond;
                  %                disp('Added *');
                  %bond
                  end 
                
                %disp('==============================');
              %end
              
          else
              error('Atom indexes must be bigger than zero.');
          end

        end
        
      end
      %%
      function add_soc(self,amp,index1,index2,varargin)
        if(self.is_soc)
          self.add_hopping(amp,index1,index2,[0 0 0],varargin{:},'mode','append');
        end
      end  
      %%
      function calculate_band(self,kvec)
        
        kx = kvec{1};
        ky = kvec{2};
        kz = kvec{3};
        
        self.E= calc_band_internal(self,kvec);
        
      end
      function Energy_cell = calc_band_internal(self,k)

        %We have NxN matrix and N determines the number of bands
        %Size of the matrix is number of atoms inside the unit cell and number of orbitals
        matrix_row_size = size(self.unit_cell,2) * self.no_orbital; % if there is soc we have 2 times bigger matrix
        if(self.is_soc), matrix_row_size = matrix_row_size * 2; end


        Energy_cell = cell(1,matrix_row_size);
        for i = 1:matrix_row_size
          Energy_cell{i} = zeros(size(k,1),1);
          %zeros(%zeros(alen,blen,clen);
        end
       
        bonding_no = numel(self.bonds);
        if(size(self.unit_cell,2) == 0)
          error('Unit cell is undefined !');
        end
        tic;
        disp('Starting diagonalization...');

        for kiter = 1:size(k,1)
          eig_matrix = zeros(matrix_row_size,matrix_row_size);

          %k = [kx(a,b,c),ky(a,b,c),kz(a,b,c)];
          for ob_iter = 1:bonding_no
              bond = self.bonds{ob_iter};
              row = bond.i;
              col = bond.j;
              kp = k(kiter,:);
              for iter_amp = 1:numel(bond.amp)
                R = bond.phase(iter_amp,:)*self.pvec;
                kdotr = 1i*(kp(1)*R(1) + kp(2)*R(2) + kp(3)*R(3));
                q = exp(kdotr);
                eig_matrix(row,col) = eig_matrix(row,col) +  (bond.amp(iter_amp,1) * q);
              end
          end
          eigens = sort(real(eig(eig_matrix)));
          for iter_eig = 1:size(eigens,1)
            Energy_cell{iter_eig}(kiter) = (eigens(iter_eig));
          end
         end

        toc;
      end

      function Energy_cell = calc_band_internalHR(self,k)

        matrix_row_size = size(self.unit_cell,2) * self.no_orbital; % if there is soc we have 2 times bigger matrix
        if(self.is_soc), matrix_row_size = matrix_row_size * 2; end


        Energy_cell = cell(1,matrix_row_size);
        for i = 1:matrix_row_size
          Energy_cell{i} = zeros(size(k,1),1);
          %Energy_cell{i} = zeros(size(kx,1),size(kx,2),size(kx,3));
        end
        
        tic;
        disp('Starting diagonalization..');
        fprintf("K iteration number :%d\n",size(k,1));

        R = self.bonds.R;
        matrix = self.bonds.matrix;
        matrix_row_size2 = matrix_row_size*matrix_row_size;
        Zeros = zeros(matrix_row_size);
        z = 1i;
        for kiter = 1:size(k,1)
          eig_matrix = Zeros;
          idpt = 1;
          kp = k(kiter,:);
          for iter = 1:size(R,1)
            Rp = R(iter,:);%*self.pvec;
            kdotr = z * (kp(1)*Rp(1) + kp(2)*Rp(2) + kp(3)*Rp(3)); 
            q = exp(kdotr);
            eig_matrix(:,:) = eig_matrix(:,:) + matrix(:,:,iter) .* q; 
          end
          eigens = sort(real(eig(eig_matrix)));
          for iter_eig = 1:size(eigens,1)
            %For some cases there are imaginary values around 1e-20 level so we will cut them with real
            Energy_cell{iter_eig}(kiter) = eigens(iter_eig) - self.Efermi;
          end
          if(mod(kiter,1000) == 0)
            fprintf("K iteration:[%d]\n",kiter);
          end
        end

      toc;
      end %function end
      %%
      function kvec = set_kvector(self,from,to,len)
        k = linspace(from/self.spatial_unit,to/self.spatial_unit,len);
        if(self.no_primvec == 3)
            [kx,ky,kz] = meshgrid(k);
        elseif(self.no_primvec == 2)
            [kx,ky] = meshgrid(k);
            kz = zeros(size(kx,1),size(kx,2));
        elseif(self.no_primvec == 1)
            [kx] = meshgrid(k);
            ky = zeros(size(kx,1),size(kx,2));
            kz = ky; 
        end    

        kk = [kx(:),ky(:),kz(:)];

        % for a = 1:size(kx,1)
        %   for b = 1:size(kx,2)
        %     for c = 1:size(kx,3)
        %       kk(end+1,:) = [kx(a,b,c),ky(a,b,c),kz(a,b,c)];
        %     end
        %   end
        % end

        kvec = {kx,ky,kz,kk};


      end
      %%
      function surfaces = plot_energy_band(self,fig,kvec,plot_type,varargin)
        %error('Dont call me we need to change k vector issue')
        %
        kx = kvec{1};
        ky = kvec{2};
        kz = kvec{3};
        k = kvec{4};

        if(isempty(self.hr_file))
          E = calc_band_internal(self,k);
        else
          E = calc_band_internalHR(self,k);
        end

        self.E = E;
        %ax = fig.CurrentAxes;
        f = fig();
        surfaces = {};
        
        kx = kx*1e-9;
        ky = ky*1e-9;
        kz = kz*1e-9;
    
        if(plot_type == 'surface')  
          hold on;
          for i = 1:size(E,2)
            if(self.no_primvec == 3) 
              Eband_sz = int32((numel(E{i}))^(1/3));
              Esurf = reshape(E{i},Eband_sz,Eband_sz,Eband_sz);
              surfaces{end+1} = surf(kx(:,:,1),ky(:,:,1),Esurf(:,:,1),varargin{:});
            elseif(self.no_primvec == 2)
              Eband_sz = int32((numel(E{i}))^(1/2));
              Esurf = reshape(E{i},Eband_sz,Eband_sz);
              surfaces{end+1} = surf(kx,ky,Esurf,varargin{:});
            elseif(self.no_primvec == 1)
              error('1d surf not implemented yet ')
                surfaces{end+1}= surf(kx,E{i},varargin{:});
            else
              error('Before plotting, define primitive vectors');
            end
            f.UserData = surfaces;
          end
        elseif(plot_type == 'contour')
          hold on;
          for i = 1:size(E,2)
            if(self.no_primvec == 3)
              error('no contour plot for 3dim ')
              %s.positive_surface = surf(kx(:,:,1),ky(:,:,1),Ep(:,:,1),varargin{:});  
            elseif(self.no_primvec == 2)
              surfaces{end+1} = contour(kx,ky,E{i},varargin{:});
            elseif(self.no_primvec == 1)
              surfaces{end+1}= contour(kx,E{i},varargin{:});
            else
              error(message('Before plotting, define primitive vectors'));
            end
          end
          f.UserData = surfaces;          
        end
        ax = f.CurrentAxes;
        xlabel(ax,sprintf('$k_{x}(%s^{-1})$',self.spatial_unit_str),'Interpreter','Latex','FontSize',15);
        ylabel(ax,sprintf('$k_{y}(%s^{-1})$',self.spatial_unit_str),'Interpreter','Latex','FontSize',15);
        zlabel(ax,'$E(eV)$','Interpreter','Latex','FontSize',15);
        set(ax,'TickLabelInterpreter','Latex');
        %title(ax,'Energy Band','Interpreter','Latex','FontSize',15);
        hold off;
      end
      %%
      function fs = plot_fermi_surface(self,fig,kvec,band_no,varargin)
        f = fig();
        kx = kvec{1};
        ky = kvec{2};
        kz = kvec{3};
        k = kvec{4};

        if(iscell(self.E) && isempty(self.E))
          if(isempty(self.hr_file))
            E = calc_band_internal(self,k);
          else
            E = calc_band_internalHR(self,k);
          end 
          self.E = E;         
        end


        E = self.E{band_no};  

        if(self.no_primvec == 3) 
          Eband_sz = int32((numel(E))^(1/3));
          Esurf = reshape(E,Eband_sz,Eband_sz,Eband_sz);
        elseif(self.no_primvec == 2)
          Eband_sz = int32((numel(E))^(1/2));
          Esurf = reshape(E,Eband_sz,Eband_sz);
        elseif(self.no_primvec == 1)
          error('1d surf not implemented yet ')
        else
          error('Before plotting, define primitive vectors');
        end

        isosurface(kx,ky,kz,Esurf,self.Efermi);
        colormap(gray);
        % if(E == 0) 
        %   error('First solve for energy eigens');
        % end

        fs  = 0;
      end      
      %%
      function plots = plot_high_symmetry_points(self,fig,precision,varargin)

        plots = {};
        
        nvargin = nargin - 3;

        if(nvargin < 2)
            error('In order to plot at least two high symmetry points are required.');
            return;
        end

        if(size(precision,2) == 1)
          precision = precision*ones(1,nvargin-1);
        end

        label_cell = {};
        if(iscell(varargin{2}))
          for i = 1:nvargin
            label_cell{end+1} = varargin{i}{2};
          end
        end

        kx = [];
        ky = [];
        kz = [];
        for i = 1:nvargin-1
            if(iscell(varargin{2}))
              point = varargin{i}{1};
              point2 = varargin{i+1}{1};
            else
              point  = varargin{i};
              point2 = varargin{i+1};
            end

            kx = [kx,linspace(point(1),point2(1),precision(i))];
            ky = [ky,linspace(point(2),point2(2),precision(i))];
            kz = [kz,linspace(point(3),point2(3),precision(i))];
        end

        k =[kx(:),ky(:),kz(:)]/self.spatial_unit;


        if(isempty(self.hr_file))
          E = calc_band_internal(self,k);
        else
          E = calc_band_internalHR(self,k);
        end

        hold on;
        for i = 1:size(E,2)
          plots{end+1} = plot(E{i},'-');
        end
        ax = fig.CurrentAxes;

        %title(ax,'High Symmetry Points','Interpreter','Latex','FontSize',15);
        
        ylabel(ax,'$$E(eV)$$','Interpreter','Latex','FontSize',15);
        xlabel(ax,'High Symmetry Points','Interpreter','Latex','FontSize',15);

        axtick = [0];
        for i = 1:numel(precision)
          axtick(end+1) = axtick(end) + precision(i); 
        end
         
        xticks(ax,axtick);
        if(isempty(label_cell) == 0), xticklabels(ax,label_cell); end
        set(ax,'TickLabelInterpreter','Latex');
        grid(ax);
        
        hold off;
        fig.UserData = plots;
        %self.Epath = E;
      end
      %%
      function ddd = plot_dos(self,fig,kvec,varargin)
        k = kvec{4};

        hold on;
        
        if(isempty(self.E))
          if(isempty(self.hr_file) )
            self.E = calc_band_internal(self,k);
          else
            self.E = calc_band_internalHR(self,k);
          end
        end

        E = self.E;

        omin = min(min(E{1}));
        omax = max(max(E{size(E,2)}));
        onum = 10*size(kvec{1},1);
        omega = linspace(omin,omax,onum);
        %size(E,2)

        sigma = (omax- omin)/onum*6;

        dosE = zeros(1,onum);
        for oiter = 1:numel(omega)
          omg = omega(oiter);
          sum_ = 0;
          for band_iter = 1:size(E,2)
            Eband = E{band_iter};
            for kiter = 1:size(k,1)
                sum_ = sum_ + delta_gaussian(omg - Eband(kiter), sigma);
            end
          end
            dosE(oiter) = dosE(oiter) + sum_;
        end

        ddd = dosE;
        ax = fig.CurrentAxes;
        plot(ax,omega(1,:),dosE(1,:),varargin{:});
        xlabel(ax,'E(eV)','Interpreter','Latex');
        ylabel(ax,'D(E)','Interpreter','Latex');
        %title(ax,'Density of States','Interpreter','Latex');
        set(ax,'TickLabelInterpreter','Latex');
        grid(ax);

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

        for i = 1:numel(self.bonds)
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

        line_object = "None";
        atom_object = "None";

        %gp.fig.CurrentAxes.XLim
        %[-self.canvas_x/2 self.canvas_x/2];

        lattice_fill_factorx = gp.fig.CurrentAxes.XLim(1): gp.fig.CurrentAxes.XLim(2);
        lattice_fill_factory = gp.fig.CurrentAxes.YLim(1): gp.fig.CurrentAxes.YLim(2);
        lattice_fill_factorz = gp.fig.CurrentAxes.ZLim(1): gp.fig.CurrentAxes.ZLim(2);
        %gp.xaxis_symmetric();
        %gp.yaxis_symmetric();
        %gp.zaxis_symmetric();

        plot_atoms = 1;
        plot_bonds = 1;
        plot_unit_vector = 0;

        constraint_x = lattice_fill_factorx;
        constraint_y = lattice_fill_factory;
        constraint_z = lattice_fill_factorz;
        
    
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
              %lattice_fill_factorx = varargin{i};
              constraint_x = varargin{i};
              %gp.change_axis_limits("x",lattice_fill_factorx(1),lattice_fill_factorx(end));
            case 'y'
              %lattice_fill_factory = varargin{i};
              constraint_y = varargin{i};
              %gp.change_axis_limits("y",lattice_fill_factory(1),lattice_fill_factory(end));
            case 'z'
              %lattice_fill_factorz = varargin{i}; 
              constraint_z = varargin{i};
              %gp.change_axis_limits("z",lattice_fill_factorz(1),lattice_fill_factorz(end));
            case 'only atoms'
              plot_bonds = 0; 
              i = i -1;
            case 'only bonds'
              plot_atoms = 0;
              i = i -1;
            case 'unit vector'
              plot_unit_vector = 1;
              uv_cx = varargin{i}(1);
              uv_cy = varargin{i}(2);
              uv_cz = varargin{i}(3);
            otherwise
              error('Undefined entry valid entries [atoms,bonds,x,y,z]');
          end
          i = i+1;
        end

        if(isstring(line_object) && strcmp(line_object,"None"))
          line_object = gp.draw("line black",0,0,0,0,'Visible','off');
        end
        if(isstring(atom_object) && strcmp(atom_object,"None"))
          atom_object = cell(1,size(self.unit_cell,2));
          colors = {"red","blue","green","magenta","yellow","black"};
          for i = 1:size(self.unit_cell,2)
            atom_object{i} = gp.draw("point "+colors{mod(i,6)+1},0,0,0.2,'Visible','off');
          end
        end


        %Bununla ilgili bir sorun var, lattice factoru  = 1 yapınca tek boyutlu zinciri cizmiyor 
        %if(self.dim == 1)
        %  lattice_fill_factory = 1;
        %end

        
        normalize = 1/self.spatial_unit;
        
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
        cz = 0;

        if(plot_atoms)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %first get the atoms from the unitcell
        uc_size = numel(self.unit_cell);

        old_x = 0;
        old_y = 0;
        old_z = 0;

        for i = 1:uc_size
          atom = self.unit_cell{i};

          for e1 = lattice_fill_factorx
            for e2 = lattice_fill_factory
              for e3 = lattice_fill_factorz

                txyz = e1.*a1 + e2.*a2 + e3.*a3;

                x = txyz(1) + atom.pos(1)*normalize;
                if(any(x < constraint_x(1)) || any(x > constraint_x(end))), continue, end

                y = txyz(2) + atom.pos(2)*normalize;
                if( any(y < constraint_y(1)) || any(y > constraint_y(end))), continue, end
                z = 0;

                if(self.dim > 2)
                  z = txyz(3) + atom.pos(3) * normalize;
                end
                
                if(any(z < constraint_z(1)) || any(z > constraint_z(end)) ),continue, end

                if(x == old_x & y == old_y & z == old_z) , continue, end

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
        end%end of plot atoms        
        plotted_bonds = {};
        
        if(plot_bonds)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BOND plottign
        old_x = 0;
        old_y = 0;
        old_z = 0;
        
        ctr = 0;
        plotted = 0;
       
        
        for i = 1 : numel(self.bonds)
        bond = self.bonds{i};
        pass = 0;
        if(bond.atoms{1} == bond.atoms{2})
            continue;
        end
        %
        for qqq = 1:numel(plotted_bonds)
            bd = plotted_bonds{qqq};
            disp('');
            if(all(all(bd.phase == bond.phase)) && bd.atoms{1} == bond.atoms{1} && bd.atoms{2} == bond.atoms{2})
                pass=1;
                break;
            end
            %check complex conjugate
            if(all(all(bd.phase == -bond.phase)) && bd.atoms{1} == bond.atoms{2} && bd.atoms{2} == bond.atoms{1}) 
                pass = 1;
                break;
            end
        end
        %
        if(pass),continue,end;
        
        for piter = 1:size(bond.phase,1) 
        %size(self.bonds,2)
        %size(phases,2)
        bphase = bond.phase(piter,:);
        
        
        if(bond.amp(piter) == 0)
            continue;
        end
        %if(all(bphase == 0) )
        %    continue;
        %end

          % tvec translates our unitcell to new unitcells
          tvector = a1 .* bphase(1);
          if(size(bphase,2) > 1)
            tvector = tvector + a2 .* bphase(2);
            if(size(bphase,2) > 2)
              tvector = tvector + a3 .* bphase(3);
            end
          end 
          
          atom1 = self.unit_cell{bond.atoms{1}};
          atom2 = self.unit_cell{bond.atoms{2}};

          for e1 = lattice_fill_factorx
          for e2 = lattice_fill_factory
          for e3 = lattice_fill_factorz

            txyz = e1.*a1 + e2.*a2 + e3.*a3;

            x = [atom1.pos(1)*normalize, tvector(1) + atom2.pos(1)*normalize] + txyz(1);
            if(any(x < constraint_x(1)) || any(x > constraint_x(end))), continue, end

            y = [atom1.pos(2) * normalize, tvector(2) + atom2.pos(2) * normalize]+ txyz(2);
            if( any(y < constraint_y(1)) || any(y > constraint_y(end))), continue, end

            if(self.dim > 2)
              z = [atom1.pos(3) * normalize, tvector(3) + atom2.pos(3) * normalize]+txyz(3);
            else 
              z = [0 0];
            end
            if(any(z < constraint_z(1)) || any(z > constraint_z(end)) ),continue, end 

            if(all(x == old_x) && all(y == old_y) && all(z == old_z))
              continue;
            end

            if(x(1) == x(2) && y(1) == y(2) &&  z(1) == z(2))
              continue;
            end
            
            gp.copy_to(line_object,x,y,z,'Visible','on');
            %gp.draw('text',mean(x),mean(y),mean(z),sprintf("%d",ctr),'FontSize',10);
            %ctr = ctr+1;
            plotted = 1;
            old_x = x;
            old_y = y;
            old_z = z;
          
          end
          end
          end
        end
          if(plotted)
            plotted_bonds{end+1} = bond;
            plotted = 0;
          end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end

        unit_str = self.spatial_unit_str;
        xs = sprintf("$$X(%s)$$",unit_str);
        ys = sprintf("$$Y(%s)$$",unit_str);
        zs = sprintf("$$Z(%s)$$",unit_str);
        grid(gp.fig.CurrentAxes);
        gp.fig.CurrentAxes.TickLabelInterpreter = 'latex';
        gp.set_xlabel(xs,'Interpreter','Latex');
        gp.set_ylabel(ys,'Interpreter','Latex');
        gp.set_zlabel(zs,'Interpreter','Latex');
        
        if(plot_unit_vector)
            v1 = gp.draw('vector black',uv_cx,uv_cy,uv_cx+a1(1),uv_cy+a1(2),'LineWidth',2,'MaxHeadSize',0.5);
            gp.draw('text',uv_cx+a1(1),uv_cy+a1(2),"$$a_{1}$$",'Interpreter','Latex','FontSize',15);
            v2 = gp.draw('vector black',uv_cx,uv_cy,uv_cx+a2(1),uv_cy+a2(2),'LineWidth',2,'MaxHeadSize',0.5);
            gp.draw('text',uv_cx+a2(1),uv_cy+a2(2),"$$a_{2}$$",'Interpreter','Latex','FontSize',15);
        end
        

      end %plot_lattice end
      %%
      function str_hamiltonian = symbolic_hamiltonian(self,varargin)
        %Create hamiltonian matrix

        %if use_exact == 1, we will use numerical values of amplitudes
        if(find(strcmp(varargin,'exact')))
          use_exact = 1;
        else
          use_exact = 0;
        end

        matrix_row_size = size(self.unit_cell,2)*self.no_orbital;
        if(self.is_soc), matrix_row_size = 2*matrix_row_size; end
        
        H = sym(zeros(matrix_row_size, matrix_row_size));

        k = sym('k'); 
        a = sym('a_%d',[1 3]);
        img = sym('i');

        bond_amp = sym('A_%d',[matrix_row_size,matrix_row_size]);

        acoef = 0;

        bonding_no = numel(self.bonds);
        for ob_iter = 1:bonding_no
          bond = self.bonds{ob_iter};
          phase = bond.phase;
          
          if(size(phase,2) == 1) 
            acoef = [phase(1) 0 0];
          elseif(size(phase,2) == 2) 
            acoef = [phase(1) phase(2) 0];
          elseif(size(phase,2) == 3)
            acoef = phase;
          end

          if(bond.symbolic_amp ~= 'None') 
            if(use_exact == 1)
              symbolic_amp = bond.amp;
            else
              symbolic_amp = bond.symbolic_amp;
            end
          else
            if(use_exact == 1)
              symbolic_amp = bond.amp;
            else
              symbolic_amp = bond_amp(bond.i,bond.j);
            end  
          end
          %disp(acoef)

          exponent = -img*k*(a(1)*acoef(1)+a(2)*acoef(2)+a(3)*acoef(3));
          exponent = simplify(exponent);

          q = exp(exponent);
          
          %if( isequal(H(bond.i,bond.j), ) )
          H(bond.i,bond.j) = H(bond.i,bond.j) + (symbolic_amp * q);
          %disp(H);
        end
        str_hamiltonian = simplify(H);


      end
      %%
      function ln = plot_brillouin_zone(self,gp,zone_no,varargin)
        %ln is line of the brillouin zone, so user can change its properties later

        if(zone_no < 1)
          error("Brillouin zone cannot be smaller than zero 1, enter valid number.");
        end

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
        
        from_to = -1:1;
        b_1 = b{1};
        b_2 = b{2}; 
        if(self.dim == 2)
          b_3 = [0,0,0];
        else
          b_3 = b{3};
        end

        if(1)
          for a1 = from_to
            for a2 = from_to
              for a3 = from_to
                p = a1.*b_1+a2.*b_2+a3.*b_3;
                points(end+1,:) = [p(1),p(2),p(3)];
                %points = [points, [p(1);p(2);p(3)] ];
              end
            end
          end

          if(0)
            %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            %Add k-space option for user so we can plot k space points
            gp.draw('point black',points(:,1),points(:,2),points(:,3));
          end

          %hp = hypot(points(1,:),points(2,:),points(3,:));
          hp  = sqrt(points(:,1).^2 + points(:,2).^2 + points(:,3).^2);
          uhp = uniquetol(hp,1e-4);

          good_index = [];

          for i = 1:size(hp,1)
            if( abs(uhp(zone_no+1)-hp(i))< 1e-4 )
              good_index(end+1) = i;
            end
          end

          good_index;

          vopt = points(good_index,:);
          vopt(end+1,:) = [0,0,0];
          vopt = unique(vopt,'row');
          zopt = vopt;

          gp.draw('point black',zopt(:,1),zopt(:,2),zopt(:,3));

          if(all(vopt(:,3) == 0))
            vopt(:,3) = [];
          end

          %useing voronoi calculate region vertices
          %old[vpx,vpy] = voronoi(points(1,:),points(2,:));
          %[vpx,vpy] =voronoi(vopt(1,:),vopt(2,:));
          %vopt
          P = voronoin(vopt);
          if(any(P(1,:) == Inf))
            P(1,:) = [];
          end

        [c,v] = voronoin(vopt);

        nx = c(v{int32(numel(v)/2+1)},:)
          if(any(nx(1,:) == Inf))
            nx(1,:) = [];
          end        
        tri = convhulln(nx)
        fh=figure(1);
        plot3(X(:,1),X(:,2),X(:,3),'b.','markersize',10);
        for i = 1:size(tri,1)
          patch(nx(tri(i,:),1),nx(tri(i,:),2),nx(tri(i,:),3),i,'FaceAlpha',0.85);
        end

        return;

          %Voronoi calculates a lot of points but chose the first (XXXTODO think about this later, is it ok ?)
          vpx1 = P(:,1)'; %vpx(1,:)
          vpy1 = P(:,2)'; %vpy(1,:)
          if(size(vopt,2) == 3)
            vpz1 = P(:,3)';
          else
            vpz1 = zeros(1,size(vpx1,2));
          end
          %P(:,1)'
          % vpx1 = vpx(1,:);
          % vpy1 = vpy(1,:);


          % %now we will eliminate same vertex points, unique is not forking as expected why?
          % i = 1;
          % while( i <= size(vpx1,2) )
          %   pt(1) = vpx1(i);
          %   pt(2) = vpy1(i);
          %   j = i+1;
          %   while(j <= size(vpx1,2) )
          %     if(pt(1) == vpx1(j) & pt(2) == vpy1(j))
          %       vpx1(j) = [];
          %       vpy1(j) = [];
          %     end
          %     j = j+1;
          %   end
          %   i = i+1;
          % end
          % %Now vpx1 and vpy1 holds different vertex points, which one is closes to the origin ?, using hypot we can find closes points

          % hypots = transpose(hypot(vpx1(:),vpy1(:)));
          % %now sort hypots and sort vpx1 and vpy1 according to hypot
          % hypots
          % vpx1
          % [hypots,I] = sort(hypots);
          % vpx1 = vpx1(I);
          % vpy1 = vpy1(I);

          % vpx1
          % hypots
          % %hypots

          % %get hypot minimum
          % min_list = (uniquetol(hypots,1e-4));
          % hypotmin = min_list(2);

          % i = 1;
          % while (i <= size(hypots,2))
          %   if(abs(hypots(i)-hypotmin) > 1e-4) % k-space is huge and truncation will cause some equality errors,our tolerance is 1e-4, 
          %     %if there is a difference between hypot values then we will eliminate them until finding the closest points
          %     hypots(i) = [];
          %     vpx1(i) = [];
          %     vpy1(i) = [];
          %   else
          %     i = i + 1;
          %   end 
          % end

          % %hypots


          % %now vpx1 and vpy1 hold closes points but in order to draw the brillouin zone we need to sort them in terms of angle between origin
          angles = atan2d(vpy1,vpx1);
          [angles,I] = sort(angles);
          vpx1 = vpx1(I);
          vpy1 = vpy1(I);

            
          if(plot_points)
            for i = 1:size(vpx1,2)
              gp.copy_to(point_obj,vpx1(i),vpy1(i),vpz1(i),'Visible','on');
            end
          end

          if(plot_coordinates)
            for i = 1:size(vpx1,2)
              x = fix(vpx1(i));
              y = fix(vpy1(i));
              z = fix(vpz1(i));
              txt = pretty_print_scientific(x,2) + "," + pretty_print_scientific(y,2);
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
            vpz1(end+1) = vpz1(1);
            %and now plot as a line
            ln = gp.copy_to(line_obj,vpx1,vpy1,vpz1,'Visible','on');
          end

        end

      end
   end %methods end
end %class end


function res = delta_gaussian(x,sigma)
  A = 1/sqrt(2*pi)/sigma;
  pow = (x.^2)./(2*sigma*sigma);
  res = A*exp(-pow);
end