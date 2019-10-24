classdef lattice_drawer < handle
   properties
       fig = 0;
       canvas_x = 0;
       canvas_y = 0;
   end
   methods
       
       function f = lattice_drawer(fig,varargin)
            f.fig = fig;
            if(nargin>1) 
                f.canvas_x = varargin{1};
                f.canvas_y = varargin{2};
                axis([0 f.canvas_x 0 f.canvas_y]);
            else
                f.canvas_x = fig.CurrentAxes.XLim(2)-fig.CurrentAxes.XLim(1);
                f.canvas_y = fig.CurrentAxes.YLim(2)-fig.CurrentAxes.YLim(1);
            end
       end
       
       function r = draw(self,type_str,varargin)
           self.fig();
           if(nargin < 3)
               disp('Draw function takes bigger than 3 arguments'); 
               return;
           end
           
           splitted_str = strsplit(type_str,' ');
           type = splitted_str{1};
           color = 0;
           if(size(splitted_str,2) > 1 )
                color = splitted_str{2};
           end
           
           x = varargin{1};
           y = varargin{2};
           w = varargin{3};
           
           if(x > self.canvas_x || y > self.canvas_y)
                r = 0;
                return;
           end
           
           switch type
                case 'rect'
                    h = varargin{4};
                    r = rectangle('Position',[x,y,w,h]);
                    if(color ~= 0) r.FaceColor = color; end
                case 'square'
                    r = rectangle('Position',[x,y,w,w]);
                    if(color ~= 0) r.FaceColor = color; end
                case 'circle'
                    %shift the circle to make x and y as center coordinates
                    %default is edge coordinate of rectangle
                    %x-r : edge x
                    x = x-w;
                    y = y-w;
                    %y-r : edge y
                    r = rectangle('Position',[x,y,2*w,2*w],'Curvature',[1 1]);
                    if(color ~= 0) r.FaceColor = color; end
               case 'line'
                    %x1,y1,x2,y2
                    y2 = varargin{4};
                    r = line([x,w],[y,y2]);
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
                   a = w;
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
                    y2 = varargin{4};
                    x2 = w;
                    
                    ax = self.fig.CurrentAxes;
                    hold on;
                    r = quiver(ax,x,y,x2-x,y2-y);
                    r.AutoScaleFactor = 1;
                    r.MaxHeadSize = 5/sqrt((x2-x)^2+(y2-y)^2);
                    
                   if(color ~= 0)
                        r.Color = color;
                   end
                    
                otherwise
                    disp('Draw function undefined type!');
                    r = 0;
                    return;
           end
           
           r.Tag = type;
           
               
                  

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
           if(a.Tag ~= 'vector')
            p = a.Position;
           end
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
               case 'vector'
                   cx = a.XData + a.UData/2;
                   cy = a.YData + a.VData/2;
                   
                    %XData: 30
                    %YData: 30
                    %ZData: []
                    %UData: 8
                    %VData: 3    
               otherwise
                   disp('Undefined tag for get_center function');
                   return;
           end
       end
       
       function save_pdf(self,pdfname)
            saveas(self.fig,pdfname+".pdf");
       end
       
       function r = set_title(self,text,varargin)
            r = title(self.fig.CurrentAxes,text,varargin{:});
       end
       function r = set_xlabel(self,text,varargin)
            r = xlabel(self.fig.CurrentAxes,text,varargin{:});
       end
       function r = set_ylabel(self,text,varargin)
            r = ylabel(self.fig.CurrentAxes,text,varargin{:});
       end       
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