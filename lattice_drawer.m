classdef lattice_drawer < handle
   properties
       fig = 0;
       canvas_x = 0;
       canvas_y = 0;
   end
   methods
       
       function f = lattice_drawer(fig_no,canvas_x,canvas_y)
            f.fig = figure(fig_no);
            axis([0 canvas_x 0 canvas_y]);
            f.canvas_x = canvas_x;
            f.canvas_y = canvas_y;
            
       end
       
       function r = draw(self,type,varargin)
           self.fig();
           if(nargin < 3)
               disp('Draw function takes bigger than 3 arguments'); 
               return;
           end
           x = varargin{1};
           y = varargin{2};
           w = varargin{3};
           fcolor_index = 0;
           ecolor_index = 0;
           switch type
                case 'rect'
                    h = varargin{4};
                    r = rectangle('Position',[x,y,w,h]);
                    if(nargin > 6) fcolor_index = 5; end
                    if(nargin > 7), ecolor_index = 6; end
                case 'square'
                    r = rectangle('Position',[x,y,w,w]);
                    if(nargin > 5) fcolor_index = 4; end
                    if(nargin > 6), ecolor_index = 5; end 
                case 'circle'
                    %shift the circle to make x and y as center coordinates
                    %default is edge coordinate of rectangle
                    %x-r : edge x
                    x = x-w;
                    y = y-w;
                    %y-r : edge y
                    r = rectangle('Position',[x,y,2*w,2*w],'Curvature',[1 1]);
                    if(nargin > 5) fcolor_index = 4; end
                    if(nargin > 6), ecolor_index = 5; end
               case 'line'
                    %x1,y1,x2,y2
                    y2 = varargin{4};
                    r = line([x,w],[y,y2]);
                    if(nargin > 6) r.Color = varargin{5}; end  
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
               case 'eqtri'
                   % x and y are the center of equilateral triangle
                   % w is one side of the triangle
                   a = w;
                   y1 = y-a/(2*sqrt(3));
                   y3 = y1;
                   y2 = y+a/sqrt(3);
                   x1 = x-a/2;
                   x3 = x+a/2;
                   x2 = x;
                   r = self.draw('tri',x1,y1,x2,y2,x3,y3);
               case 'hexagon'
                   a = w;
                   x1 = x-a; y1 = y;
                   x2 = x-a/2; y2 = y+a*sqrt(3)/2;
                   x3 = x+a/2; y3 = y2;
                   x4 = x+a; y4 = y;
                   x5 = x3; y5 = y-a*sqrt(3)/2;
                   x6 = x2; y6 = y5;
                   
                   r = self.draw('line',x1,y1,x2,y2);
                   r.UserData(end+1) = self.draw('line',x2,y2,x3,y3);
                   r.UserData(end+1) = self.draw('line',x3,y3,x4,y4);
                   r.UserData(end+1) = self.draw('line',x4,y4,x5,y5);
                   r.UserData(end+1) = self.draw('line',x5,y5,x6,y6);
                   r.UserData(end+1) = self.draw('line',x6,y6,x1,y1);
               case 'vector'
                    y2 = varargin{4};
                    r = line([x,w],[y,y2]);
                    r.Marker = 'o';
                    r.MarkerIndices = [2 2];
                    if(nargin > 6) 
                        r.Color = varargin{5};
                        r.MarkerFaceColor = r.Color;
                    end
                    
                otherwise
                    disp('Draw function undefined type!');
                    r = 0;
                    return;
           end
           
           r.Tag = type;
           if(fcolor_index ~= 0) r.FaceColor = varargin{fcolor_index}; end
           if(ecolor_index ~= 0) r.EdgeColor = varargin{ecolor_index}; end
               
                  

       end
        
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
       
       function set_text(self,obj,txt,varargin)
           [cax,cay] = self.get_center(obj);
           if(nargin>3)
            text(cax,cay,txt,varargin{:});
           else
               text(cax,cay,txt);
           end
       end
       
       function [cx,cy] = get_center(self,a)
           p = a.Position;
           switch a.Tag
               case 'rect'
                   cx = (p(1)+p(3)/2);
                   cy = (p(2)+p(4)/2);
               case 'square'
                   cx = (p(1)+p(3)/2);
                   cy = (p(2)+p(4)/2);
               case 'circle'
                   r = p(3)/2;
                   cx = p(1)+r;
                   cy = p(2)+r;
               otherwise
                   disp('Undefined tag for get_center function');
                   return;
           end
       end
       
       
   end
end

%axis([0 100 0 100])
%pos = [50 50 10 10]; 
%r = rectangle('Position',pos,'Curvature',[1 1]);
