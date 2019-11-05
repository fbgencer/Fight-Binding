classdef lattice_drawer < handle
	 properties
			 fig = 0;
			 canvas_x = 0;
			 canvas_y = 0;
			 canvas_z = 0;
	 end
	 methods
			 %%
			 function f = lattice_drawer(fig,varargin)
						f.fig = fig;
						nvargin = nargin-1;
						if(nvargin > 1) 
								f.canvas_x = varargin{1};
								f.canvas_y = varargin{2};
								if(nvargin > 2)
									f.canvas_z = varargin{3};
								end
								%fig.CurrentAxes.XLim = [0 f.canvas_x];
								%fig.CurrentAxes.YLim = [0 f.canvas_y];
								%fig.CurrentAxes.ZLim = [0 f.canvas_z];
								
								axis([0 f.canvas_x 0 f.canvas_y ]);
								if(f.canvas_z ~= 0)
									fig.CurrentAxes.ZLim = [0 f.canvas_z];
									view(fig.CurrentAxes,3);
								end

						else
								f.canvas_x = fig.CurrentAxes.XLim(2) - fig.CurrentAxes.XLim(1);
								f.canvas_y = fig.CurrentAxes.YLim(2) - fig.CurrentAxes.YLim(1);
								f.canvas_z = fig.CurrentAxes.ZLim(2) - fig.CurrentAxes.ZLim(1);
						end
			 end
			 %%
			 function change_axis_limits(self,varargin)
			 	if(nargin > 1)
			 		i = 1;
			 		while(i <= size(varargin,2))
			 			if(varargin{i} == "x")
			 				self.fig.CurrentAxes.XLim = [varargin{i+1} varargin{i+2}];
			 			elseif(varargin{i} == "y")
			 				self.fig.CurrentAxes.YLim = [varargin{i+1} varargin{i+2}];
			 			elseif(varargin{i} == "z")
			 				self.fig.CurrentAxes.ZLim = [varargin{i+1} varargin{i+2}];
			 			else
			 				disp('Coordinate strings must be x,y or z');
			 				return;
			 			end
			 			i = i+3;
			 		end
			 	end
			 end
			 %%
			 function xaxis_symmetric(self)
			 		self.fig.CurrentAxes.XLim = [-self.canvas_x/2 self.canvas_x/2];
			 end
			 %%
			 function yaxis_symmetric(self)
			 		self.fig.CurrentAxes.YLim = [-self.canvas_y/2 self.canvas_y/2];
			 end
			 %%
			 function zaxis_symmetric(self)
			 	if(self.canvas_z > 0)
			 		self.fig.CurrentAxes.ZLim = [-self.canvas_z/2 self.canvas_z/2];
			 	end
			 end			 
			 %%
			 function r = draw(self,type_string,varargin)
					 self.fig();
					 ax = self.fig.CurrentAxes;
					 if(nargin < 3)
							 disp('Draw function takes at least 2 arguments'); 
							 return;
					 end
					 
					 splitted_str = strsplit(type_string,' ');
					 type = splitted_str{1};
					 color = 0;
					 if(size(splitted_str,2) > 1 )
								color = splitted_str{2};
					 end
					 
					 x = varargin{1};
					 y = varargin{2};
					 

					 %This uncomment is important, later decide it is useful or not
					 %some points can be outside the canvas but if we return 0 text or other functions
					 %raise error
					% if(x > self.canvas_x | y > self.canvas_y)
					%			r = 0;
					%			return;
					 %end

					 %number of varargin, basically size 
					 nvargin = size(varargin,2);
					 
					 switch type
					 			case 'crect'
					 				%centered rectangle
				 					w = varargin{3};
									h = varargin{4};
									r = rectangle('Position',[x-w/2,y-h/2,w,h],varargin{5:end});
									if(color ~= 0) r.FaceColor = color; end					 				
								case 'rect'
					 					w = varargin{3};
										h = varargin{4};
										r = rectangle('Position',[x,y,w,h],varargin{5:end});
										if(color ~= 0) r.FaceColor = color; end
								case 'square'
					 					w = varargin{3};
										r = rectangle('Position',[x,y,w,w],varargin{5:end});
										if(color ~= 0) r.FaceColor = color; end
								case 'circle'

					 					rad = varargin{3};
										%shift the circle to make x and y as center coordinates
										%default is edge coordinate of rectangle
										%x-r : edge x
										x = x-rad;
										y = y-rad;
										%y-r : edge y
										r = rectangle('Position',[x,y,2*rad,2*rad],'Curvature',[1 1],varargin{4:end});
										if(color ~= 0) r.FaceColor = color; end
							 case 'line'
										
							 			if( isscalar(varargin{1}) == 0 & isscalar(varargin{2}) == 0 )
											x = varargin{1};
											y = varargin{2};
											if(isscalar(varargin{3}) == 0)
												z = varargin{3};
												varargin_start = 4;
											else
												varargin_start = 3;
												z = zeros(1,size(x,2));
											end
										elseif(nvargin > 4 & isnumeric(varargin{5}))
												%x1,y1,z1,x2,y2,z3;
												x = [ varargin{1}, varargin{4} ];
												y = [ varargin{2}, varargin{5} ];
												z = [ varargin{3}, varargin{6} ];
												varargin_start = 7;
										else
											%x1,y1,x2,y2
											x = [ varargin{1}, varargin{3} ];
											y = [ varargin{2}, varargin{4} ];
											z = [ 0,0 ];
											varargin_start = 5;
										end

										r = line(ax,x,y,z,varargin{varargin_start:end});

										if(color ~= 0) r.Color = color; end
							 case 'tri'
									 if(nargin ~= 8)
											 error(message('Triangle requires 6 coordinates'));
												return; 
									 end
									 %x1 y1 x2 y2 x3 y3
									 r = line([varargin{1} varargin{3}],[varargin{2} varargin{4} ]); % x1 x2 y1 y2
									 r2 = line([varargin{3} varargin{5}],[varargin{4} varargin{6}]); % x2 x3 y2 y3
									 r3 = line([varargin{5} varargin{1}],[varargin{6} varargin{2}]); % x3 x1 y3 y1
									 r.UserData = [r2 r3]; % we must get vertices from r.User Data
									 if(color ~= 0) 
											 r.Color = color; 
											 r2.Color = color;
											 c3.Color = color;
									 end
							 case 'eqtri'
									 % x and y are the center of equilateral triangle
									 % w is one side of the triangle
									 w = varargin{3};
									 a = w;
									 y1 = y-a/(2*sqrt(3));
									 y3 = y1;
									 y2 = y+a/sqrt(3);
									 x1 = x-a/2;
									 x3 = x+a/2;
									 x2 = x;
									 if(color ~= 0)
												r = self.draw('tri'+' '+color,x1,y1,x2,y2,x3,y3);
									 else
											 r = self.draw('tri',x1,y1,x2,y2,x3,y3);
									 end
							 case 'hexagon'
									 a = varargin{3};
									 x1 = x-a; y1 = y;
									 x2 = x-a/2; y2 = y+a*sqrt(3)/2;
									 x3 = x+a/2; y3 = y2;
									 x4 = x+a; y4 = y;
									 x5 = x3; y5 = y-a*sqrt(3)/2;
									 x6 = x2; y6 = y5;
									 
									 tp = 'line';
									 if(color ~= 0 ) tp = tp + ' ' + color; end
									 r = self.draw(tp,x1,y1,x2,y2);
									 n1 = self.draw(tp,x2,y2,x3,y3);
									 n2 = self.draw(tp,x3,y3,x4,y4);
									 n3 = self.draw(tp,x4,y4,x5,y5);
									 n4 = self.draw(tp,x5,y5,x6,y6);
									 n5 = self.draw(tp,x6,y6,x1,y1);
									 r.UserData = [r n1 n2 n3 n4 n5];

							 case 'vector'
									x2 = varargin{3};
									y2 = varargin{4};

									hold on;

									if(nvargin > 4)
										if( isnumeric(varargin{5}) )
											x = varargin{1};
											y = varargin{2};
											z = varargin{3};
											x2 = varargin{4};
											y2 = varargin{5};
											z2 = varargin{6};
											r = quiver3(ax,x,y,z,x2-x,y2-y,z2-z,varargin{7:end});
										else
											r = quiver(ax,x,y,x2-x,y2-y,varargin{5:end});
										end
									else
											r = quiver(ax,x,y,x2-x,y2-y,varargin{5:end});
									end

									r.AutoScaleFactor = 1;
									r.MaxHeadSize = 5/sqrt((x2-x)^2+(y2-y)^2);
										
									if(color ~= 0)
										r.Color = color;
									end

								case 'cuboid'
									x1 = varargin{1};
									y1 = varargin{2};
									z1 = varargin{3};
									w = varargin{4};
									l = varargin{5};
									h = varargin{6};

									vertx = [x1, x1, x1, x1, x1+w, x1+w, x1+w, x1+w];
									verty = [y1, y1+l, y1+l, y1, y1, y1+l, y1+l, y1];
									vertz = [z1, z1, z1+h, z1+h, z1, z1, z1+h, z1+h];
									fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
									r = patch(ax,vertx,verty,vertz,'red','Faces',fac,varargin{7:end});
									if(color ~= 0)
										r.FaceColor = color;
									end
									view(ax,3);

								case 'point'
									z = 0;
									varargin_start = 3;
									if(nvargin > 2)
										if(isnumeric(varargin{3}))
											z = varargin{3};
											varargin_start = 4;
										end
									end
									hold on;
									
									r = scatter3(ax,x,y,z,varargin{varargin_start:end});

									if(color ~= 0)
										r.MarkerFaceColor = color;
										r.MarkerEdgeColor = color;
									else
										r.MarkerFaceColor = 'black';
										r.MarkerEdgeColor = 'black';
									end
								case 'text'
									x = varargin{1};
									y = varargin{2};
									txt = varargin{3};

									r = text(x,y,txt,varargin{4:end});
								otherwise
										disp('Draw function undefined type!');
										r = 0;
										return;
					 end
					 
					 r.Tag = type;
					 
			 end
				%%
			 function connect(self,a,b,varargin)
						%here a and b are two different object in the plot
						%this function will draw a line between them
						[cax,cay] = self.get_center(a);
						[cbx,cby] = self.get_center(b);
						color = 'k';
						if(nargin>3)
								color = varargin{1};
						end
						self.draw('line',cax,cay,cbx,cby,color);  
			 end
			 %%
			 function t = set_text(self,obj,txt,varargin)
					[cax,cay] = self.get_center(obj);
					if(nargin>3)
							t = text(cax,cay,txt,varargin{:});
					else
							t = text(cax,cay,txt);
					end
					t.FontSize = 7;
			 end
			 %%
			 function [cx,cy] = get_center(self,a)

			 			if(a == 0)
			 				disp("Empty value in get_center")
			 				return;
			 			end

					 switch a.Tag
							case "rect"
								p = a.Position;
								cx = (p(1)+p(3)/2);
								cy = (p(2)+p(4)/2);
							case "square"
								p = a.Position;
								cx = (p(1)+p(3)/2);
								cy = (p(2)+p(4)/2);
							case "circle"
								p = a.Position;
								r = p(3)/2;
								cx = p(1)+r;
								cy = p(2)+r;
							case "line"
								cx = 0.5*(a.XData(1) + a.XData(2));
								cy = 0.5*(a.YData(1) + a.YData(2));									
							case "vector"
								cx = a.XData + a.UData/2;
								cy = a.YData + a.VData/2;
							case 'point'
								cx = a.XData;
								cy = a.YData;
							otherwise
								disp('Undefined tag for get_center function');
								return;
					 end
			 end
			 %%
			 function r = copy_to(self,obj,varargin)
						r = copy(obj);
						x = varargin{1};
						y = varargin{2};
						w = 0;

						if(nargin > 4)
							is_entry_numeric = (isstring(varargin{3}) == 0  && ischar(varargin{3}) == 0 );
						end
						nvargin = size(varargin,2);


						switch obj.Tag
								case 'crect'
									if(nargin > 4 &&  is_entry_numeric)
											w = varargin{3};
											h = varargin{4};
											set(r,varargin{5:end});
									else
											w = r.Position(3);
											h = r.Position(4);
											set(r,varargin{3:end});
									end
										r.Position = [x-w/2,y-h/2,w,h];

								case 'rect'
										if(nargin > 4 &&  is_entry_numeric)
												w = varargin{3};
												h = varargin{4};
												set(r,varargin{5:end});
										else
												w = r.Position(3);
												h = r.Position(4);
												set(r,varargin{3:end});
										end
										r.Position = [x,y,w,h];
										
								case 'circle'
										if(nargin > 4 && is_entry_numeric)
												w = varargin{3};
												set(r,varargin{4:end});
										else
												w = r.Position(3)/2;
												set(r,varargin{3:end});
										end
										x = x-w;
										y = y-w;
										r.Position = [x,y,2*w,2*w];
										
								case 'line'

									varargin_start = 1;
						 			if( isscalar(varargin{1}) == 0 & isscalar(varargin{2}) == 0 )
										x = varargin{1};
										y = varargin{2};
										if(isscalar(varargin{3}) == 0 & ischar(varargin{3}) == 0 & isstring(varargin{3}) == 0)
											z = varargin{3};
											varargin_start = 4;
										else
											varargin_start = 3;
											z = r.ZData(1)*ones(1,size(x,2));
										end
									elseif(nvargin > 4 & isnumeric(varargin{5}))
											%x1,y1,z1,x2,y2,z3;
											x = [ varargin{1}, varargin{4} ];
											y = [ varargin{2}, varargin{5} ];
											z = [ varargin{3}, varargin{6} ];
											varargin_start = 7;
									else
										%x1,y1,x2,y2
										x = [ varargin{1}, varargin{3} ];
										y = [ varargin{2}, varargin{4} ];
										z = [r.ZData(1), r.ZData(2)];
										varargin_start = 5;
									end

									r.XData = x;
									r.YData = y;
									r.ZData = z;
									% if(nvargin > 4)
									% 	if( isnumeric(varargin{5}) )
									% 		r.XData = [varargin{1},varargin{4}];
									% 		r.YData = [varargin{2},varargin{5}];
									% 		r.ZData = [varargin{3},varargin{6}];
											
									% 	else
									% 		set(r,varargin{5:end});
									% 	end
									% end
									if(nvargin >= varargin_start)
										set(r,varargin{varargin_start:end});
									end
								case 'point'
									r.XData = varargin{1};
									r.YData = varargin{2};

									varargin_start = 3;
									if(nvargin >= 3 & ischar(varargin{3}) == 0 & isstring(varargin{3}) == 0)
										r.ZData = varargin{3};
										varargin_start = 4;
									end

									if(nvargin >= varargin_start)
										set(r,varargin{varargin_start:end});
									end
						end

						r.Parent = obj.Parent;
			 end
			 %%
			 function save_pdf(self,pdfname)
						saveas(self.fig,pdfname+".pdf");
			 end
			 %%
			 function r = set_title(self,text,varargin)
						r = title(self.fig.CurrentAxes,text,varargin{:});
			 end
			 %%
			 function r = set_xlabel(self,text,varargin)
						r = xlabel(self.fig.CurrentAxes,text,varargin{:});
			 end
			 %%
			 function r = set_ylabel(self,text,varargin)
						r = ylabel(self.fig.CurrentAxes,text,varargin{:});
			 end
			 %%
			 function [rx,ry] = rotate(self,x,y,angle)
				t  = angle*pi/180;
				Rx = [cos(t) -sin(t); sin(t) cos(t)];
				a = Rx*[x;y];
				rx = a(1,1);
				ry = a(2,1);
			 end
			 
	 end
end

%axis([0 100 0 100])
%pos = [50 50 10 10]; 
%r = rectangle('Position',pos,'Curvature',[1 1]);
