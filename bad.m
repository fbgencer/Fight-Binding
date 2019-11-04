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
                atom = self.unit_cell{i}; %unit_cell{} contains atoms inside the unit cell
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
                tvector = a1 .* bond.phase(1) + a2 .* bond.phase(2); % tvec translates our unitcell to new unitcells
                new_unitcell_cx = cx + tvector(1);
                new_unitcell_cy = cy + tvector(2);
                %gp.draw('vector',cx,cy,new_unitcell_cx,new_unitcell_cy);
                x(1) = cx + atom1.pos(1) * normalize;
                y(1) = cy + atom1.pos(2) * normalize;
                x(2) = new_unitcell_cx + atom2.pos(1) * normalize;
                y(2) = new_unitcell_cy + atom2.pos(2) * normalize;
                
                if(x(1) == x(2) && y(1) == y(2))
                   continue;
                end

                tbond = gp.copy_to(plot_bond_obj{1},x(1),y(1),x(2),y(2),'Visible','on');
                for site_iter = 1:size(plot_atoms_obj,2)
                  gp.copy_to(plot_atoms_obj{site_iter},x(site_iter),y(site_iter),'Visible','on');
                end
            end
            end
            end
      end