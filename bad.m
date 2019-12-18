%Based on QHULL DEMO
%http://www.mathworks.com/products/demos/shipping/matlab/qhulldemo.html

clear all;
close all;
%clc;

disp('+++++++++++++++++++');
disp('  BRILLOUIN ZONES  ');
disp('(c)  Poms 2009     ');
disp('+++++++++++++++++++');

% a=input(['\n','General triclinic lattice','\n','a: ']);
% if isempty(a), a=1; end
% b=input('b: ');
% if isempty(b), b=1; end
% c=input('c: ');
% if isempty(c), c=1; end
% alpha=input('alpha: ');
% if isempty(alpha), alpha=72; end
% beta=input('beta: ');
% if isempty(beta), beta=80; end
% gamma=input('gamma: ');
% if isempty(gamma), gamma=110; end
% animation=input('Make an animation and save the movie as gif yes/no: ','s');
% if isempty(animation), animation='n'; end

animation = 'n';
alpha = 0;
beta = 0;
gamma = 0;

alpha=alpha*pi/180;
beta=beta*pi/180;
gamma=gamma*pi/180;

w1 = cos(alpha) - cos(beta) * cos(gamma)/(sin(beta) * sin(gamma));
w2 = sin(gamma)^2 - cos(beta)^2 - cos(alpha)^2 ...
    + 2*cos(alpha) *cos(beta) * cos(gamma);
w2 = sqrt(w2) / (sin(beta) * sin(gamma));

tp = 'fcc'


%T1 = [a, 0, 0];
%T2 = [b * cos(gamma), b * sin(gamma), 0];
%T3 = [c * cos(beta), c * w1*sin(beta), c*w2*sin(beta)];

if(tp == 'cub')
a = 1;
T1 = a*[1 0 0];
T2 = a*[0 1 0];
T3 = a*[0 0 1];
elseif(tp == 'bcc')
a = 1/2;
T1 = a*[1,1,-1];
T2 = a*[-1,1,1];
T3 = a*[1,-1,1];
else
a = 1;
T1 = a*[1 1 0];
T2 = a*[0 1 1];
T3 = a*[1 0 1];
end


spat=cross(T1,T2)*T3';

G1=2*pi*cross(T2,T3)/spat;
G2=2*pi*cross(T3,T1)/spat;
G3=2*pi*cross(T1,T2)/spat;

X = [];

from_to = -2:2;
for i=from_to
    for j=from_to
        for k=from_to
            X(end+1,:)=i.*G1+j.*G2+k.*G3;
            %fprintf("[%d,%d,%d]=(%g,%g,%g)\n",i,j,k,X(end,:));
        end
    end
end

if(1)
points = X;
hp  = sqrt(points(:,1).^2 + points(:,2).^2 + points(:,3).^2);
uhp = uniquetol(hp,1e-4);

zone_no = 5;
good_index = [];

for i = 1:size(hp,1)
if( abs(uhp(zone_no+1)-hp(i))< 1e-4 )
  good_index(end+1) = i;
end
end
vopt = points(good_index,:);
%vopt(end+1,:) = [0,0,0];
vopt = unique(vopt,'row');
%X = vopt;
end

% X(7,:) = [];
% X(8,:) = [];
% X(19,:) = [];
% X(22,:) = [];
% X(3,:) = [];
XX = X;
for j = 1:size(vopt,1)
    i = 1;
    while(i <= size(XX,1))
        if(vopt(j,:) == XX(i,:))
            %XX(i,:)
            XX(i,:) = [];

            %fprintf("found at %d %d\n",i,j);
        else
            i = i + 1;
        end
    end
end

X = XX;
cla reset; hold on

% Compute Voronoi diagram.
[c,v] = voronoin(X);

iter = floor(numel(v)/2)+1;
nx = c(v{iter},:);
tri = convhulln(nx);
if(c(1,:) == Inf)
    c(1,:)  = [];
end
%l = line(a,b,c,'LineWidth',2,'Color','k')
%bendeki nx


p = {};
if (animation ~= 'y')
    fh=figure(1);
    %plot3(c(:,1),c(:,2),c(:,3),'g.','markersize',20)
    %kpoints
    %plot3(X(:,1),X(:,2),X(:,3),'b.','markersize',10);
    %vertec pt
    %plot3(nx(:,1),nx(:,2),nx(:,3),'r.','markersize',20)
    %
    for i = 1:size(tri,1)
        p{end+1} = patch(nx(tri(i,:),1),nx(tri(i,:),2),nx(tri(i,:),3),i,'FaceAlpha',0.2,'LineStyle','none');
    end
    % Modify the view.
    view(3);
    title('1st Brillouin zone');
    %axis equal tight off  vis3d
    grid on;
    %camzoom(2);
    rotate3d on
else
    fh=figure('Position',[10,10,400,300]);  %size of the gif!
    for i = 1:size(tri,1)
        patch(nx(tri(i,:),1),nx(tri(i,:),2),nx(tri(i,:),3),i,'FaceAlpha',0.85);
    end
    title('1st Brillouin zone');
    axis equal tight off vis3d
    camtarget([0,0,0]);
    loop=36;
    azimuth=360/loop*[0:loop-1];
    view(azimuth(1),30);
    drawnow;
    pause(0.1);
    mov=getframe(gcf);
    
    gifname = 'animation.gif';
    I = frame2im(mov);
    [X, map] = rgb2ind(I, 128);
    imwrite(X, map, gifname, 'GIF', 'WriteMode', 'overwrite', 'DelayTime', 0, 'LoopCount', Inf);
    
    for k=2:loop
     view(azimuth(k),30);
      drawnow;
      pause(0.1);
      mov=getframe(gcf);
      I = frame2im(mov);
        [X, map] = rgb2ind(I, 128);
        imwrite(X, map, gifname, 'GIF', 'WriteMode', 'append', 'DelayTime', 0);
    end
end




for iter = 1:size(nx,1)
pt = nx(iter,:);
dif = nx-pt;
hpdif = sqrt(dif(:,1).^2 + dif(:,2).^2 + dif(:,3).^2);
uhpdif = uniquetol(hpdif,1e-4); %first is zero strt from 2nd
good_index2 = [];
for i = 1:size(hpdif,1)
if( abs(uhpdif(2)-hpdif(i))< 1e-4 )
  good_index2(end+1) = i;
end
end

good_points = nx(good_index2,:);
for i = 1:size(good_points,1)
    line([pt(1),good_points(i,1)],[pt(2),good_points(i,2)],[pt(3),good_points(i,3)],'LineWidth',2);
end
end