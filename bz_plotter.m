clear all;
close all;


tp = 'nbse2'


%T1 = [a, 0, 0];
%T2 = [b * cos(gamma), b * sin(gamma), 0];
%T3 = [c * cos(beta), c * w1*sin(beta), c*w2*sin(beta)];

if(strcmp(tp,'cub'))
a = 1;
T1 = a*[1 0 0];
T2 = a*[0 1 0];
T3 = a*[0 0 1];
elseif(strcmp(tp,'bcc'))
a = 1/2;
T1 = a*[1,1,-1];
T2 = a*[-1,1,1];
T3 = a*[1,-1,1];
elseif(strcmp(tp,'fcc'))
a = 1;
T1 = a*[1 1 0];
T2 = a*[0 1 1];
T3 = a*[1 0 1];
elseif(strcmp(tp,'nbse2'))
a = 1;
T1 = a*[2.9439482089090570 -1.6996892908939607 0.0000000000000000];
T2 = a*[0.0000000000000000  3.3993785817879214   0.0000000000000000];
T3 = a*[0.0000000000000000  0.0000000000000000 24.0320930230956016];

elseif(strcmp(tp,'bi2se3'))
T1 = [-2.069  -3.583614  0.000000];
T2 = [2.069  -3.583614  0.000000];
T3 = [0.000   2.389075  9.546667];

end


spat=cross(T1,T2)*T3';

G1=2*pi*cross(T2,T3)/spat;
G2=2*pi*cross(T3,T1)/spat;
G3=2*pi*cross(T1,T2)/spat;

Kvec = [G1;G2;G3];


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
if (1)
    fh=figure(1);
    gp = lattice_drawer(fh,10,10,10);
    gp.axis_symmetric();


    gp.set_xlabel('$$K_{x}(\AA^{-1})$$','Interpreter','Latex','FontSize',20);
    gp.set_ylabel('$$K_{y}(\AA^{-1})$$','Interpreter','Latex','FontSize',20);
    gp.set_zlabel('$$K_{z}(\AA^{-1})$$','Interpreter','Latex','FontSize',20);        


    G = {[0 0 0],'$$\Gamma$$'};
M = {[0.5 0 0]*Kvec,'M'};
K = {[0.33 0.33 0]*Kvec,'K'};
L = {[0.5 0 0.5]*Kvec,'L'};
    gmp = gp.draw('point black',G{1}(1),G{1}(2),G{1}(3));
    tg=gp.set_text(gmp,'$$\Gamma$$','Interpreter','Latex');
    tg.FontSize = 16;

    mp = gp.draw('point black',M{1}(1),M{1}(2),M{1}(3));
    tm=gp.set_text(mp,'M','Interpreter','Latex');
    tm.FontSize = 16;

    kp = gp.draw('point black',K{1}(1),K{1}(2),K{1}(3));
    tk=gp.set_text(kp,'K','Interpreter','Latex');
    tk.FontSize = 16;

    %kvectors
    % for i = 1:3
    % vcec1 = gp.draw('vector black',0,0,0,Kvec(i,1),Kvec(i,2),Kvec(i,3));
    % vcec1.MaxHeadSize = 0.3;
    % vcec1.LineWidth = 2;
    % end
    % vcec1.MaxHeadSize = 0.9;


    %path vectors
    vec1 = gp.draw('vector rgb:ff0000',G{1}(1),G{1}(2),G{1}(3),M{1}(1),M{1}(2),M{1}(3));
    vec1.MaxHeadSize = 0.5;
    vec1.LineWidth = 2.2;

    vec1 = gp.draw('vector rgb:ff0000',M{1}(1),M{1}(2),M{1}(3),K{1}(1),K{1}(2),K{1}(3));
    vec1.MaxHeadSize = 0.5;
    vec1.LineWidth = 2.2;

    vec1 = gp.draw('vector rgb:ff0000',K{1}(1),K{1}(2),K{1}(3),G{1}(1),G{1}(2),G{1}(3));
    vec1.MaxHeadSize = 0.5;
    vec1.LineWidth = 2.2;



    %plot3(c(:,1),c(:,2),c(:,3),'g.','markersize',20)
    %kpoints
    %plot3(X(:,1),X(:,2),X(:,3),'b.','markersize',10);
    %vertec pt
    %plot3(nx(:,1),nx(:,2),nx(:,3),'r.','markersize',20)
    %

    xy_points = {};
    vertex_points = [];

    for i = 1:size(tri,1)
        p{end+1} = patch(nx(tri(i,:),1),nx(tri(i,:),2),nx(tri(i,:),3),i,'FaceAlpha',0.1,'LineStyle','none',...
            'FaceColor',sscanf("FF0000",'%2x%2x%2x',[1 3])/255);
        px = p{end};
        %now try to build a 2d path
        vertex_points(end+1,:) = px.Vertices(1,:);
        vertex_points(end+1,:) = px.Vertices(2,:);
        vertex_points(end+1,:) = px.Vertices(3,:);

        % i
        % px.Vertices
        % point1 = px.Vertices(1,:)
        % point2 = px.Vertices(2,:)
        % point3 = px.Vertices(3,:)
        
        % abs(point1(3)-point2(3))<1e-4
        % abs(point1(3)-point3(3))<1e-4
        % abs(point2(3)-point3(3))<1e-4
        
        % if(abs(point1(3)-point2(3))<1e-4 )
        %     v1 = point1;
        %     v2 = point2;
        % elseif(abs(point1(3)-point3(3))<1e-4 )
        %     v1 = point1;
        %     v2 = point3;
        % elseif(abs(point2(3)-point3(3))<1e-4 )
        %     v1 = point2;
        %     v2 = point3;
        % else 
        %     continue;
        % end
           
        % xy_points{end+1} = [v1;v2];

        % xx = [v1(1) v2(1)];
        % xy = [v1(2) v2(2)];
        % xz = [v1(3) v2(3)];
        % q = gp.draw('line black',xx,xy,xz);
        % gp.set_text(q,sprintf('%d',i))
        % disp("=====================0")

    end

    %burayi nbse2 için aç
    if(1)
    vertex_points = unique(vertex_points,'row');
    %sort x y
    angles = atan2d(vertex_points(:,2),vertex_points(:,1));
    [angles,I] = sort(angles);
    vertex_points = vertex_points(I,:);
    
    unz = uniquetol(vertex_points(:,3),1e-4);
    faces = {};
    for i = 1:numel(unz)
        q = [];
        dm = abs(vertex_points(:,3)-unz(i))<1e-4;
        for j = 1:numel(dm)
            if(dm(j) == 1)
                q(end+1,:) = vertex_points(j,:);       
            end
        end
        faces{end+1} = q;
    end 


    for i = 1:numel(faces)
        surf1 = faces{i};
        for j = 1:size(surf1,1)
            v1 = surf1(j,:); v2 = surf1(mod(j,size(surf1,1))+1,:);
            xx = [v1(1) v2(1)];
            xy = [v1(2) v2(2)];
            xz = [v1(3) v2(3)];
            q = gp.draw('line',xx,xy,xz,'LineWidth',2);
        end
    end

    for i = 1:size(vertex_points,1)
        gp.draw('point',vertex_points(i,1),vertex_points(i,2),vertex_points(i,3));
    end

    end


    % Modify the view.
    view(3);
    %axis equal tight off  vis3d
    grid on;
    %camzoom(2);
    rotate3d on
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



%nbse2 paths
if(0)

    G = {[0 0 0],'$$\Gamma$$'};
M = {[0.5 0 0]*Kvec,'M'};
K = {[0.33 0.33 0]*Kvec,'K'};
L = {[0.5 0 0.5]*Kvec,'L'};
    gmp = gp.draw('point black',G{1}(1),G{1}(2),G{1}(3));
    tg=gp.set_text(gmp,'$$\Gamma$$','Interpreter','Latex');
    tg.FontSize = 16;

    mp = gp.draw('point black',M{1}(1),M{1}(2),M{1}(3));
    tm=gp.set_text(mp,'M','Interpreter','Latex');
    tm.FontSize = 16;

    kp = gp.draw('point black',K{1}(1),K{1}(2),K{1}(3));
    tk=gp.set_text(kp,'K','Interpreter','Latex');
    tk.FontSize = 16;

    %kvectors
    % for i = 1:3
    % vcec1 = gp.draw('vector black',0,0,0,Kvec(i,1),Kvec(i,2),Kvec(i,3));
    % vcec1.MaxHeadSize = 0.3;
    % vcec1.LineWidth = 2;
    % end
    % vcec1.MaxHeadSize = 0.9;


    %path vectors
    vec1 = gp.draw('vector rgb:ff0000',G{1}(1),G{1}(2),G{1}(3),M{1}(1),M{1}(2),M{1}(3));
    vec1.MaxHeadSize = 0.5;
    vec1.LineWidth = 2.2;

    vec1 = gp.draw('vector rgb:ff0000',M{1}(1),M{1}(2),M{1}(3),K{1}(1),K{1}(2),K{1}(3));
    vec1.MaxHeadSize = 0.5;
    vec1.LineWidth = 2.2;

    vec1 = gp.draw('vector rgb:ff0000',K{1}(1),K{1}(2),K{1}(3),G{1}(1),G{1}(2),G{1}(3));
    vec1.MaxHeadSize = 0.5;
    vec1.LineWidth = 2.2;
end

%gaas path
if(0)
    L = {[0.5 0.5 0.5]*Kvec,'L'};
    G = {[0 0 0],'$$\Gamma$$'};
    X = {[0.5 0 0.5]*Kvec,'X'};    
    gmp = gp.draw('point black',G{1}(1),G{1}(2),G{1}(3));
    tg=gp.set_text(gmp,'$$\Gamma$$','Interpreter','Latex');
    tg.FontSize = 16;
    tg.Position = [-0.2 -0.3 -0.4]

    mp = gp.draw('point black',L{1}(1),L{1}(2),L{1}(3));
    tm=gp.set_text(mp,'L','Interpreter','Latex');
    tm.FontSize = 16;

    kp = gp.draw('point black',X{1}(1),X{1}(2),X{1}(3));
    tk=gp.set_text(kp,'X','Interpreter','Latex');
    tk.FontSize = 16;


    vec1 = gp.draw('vector rgb:ff0000',L{1}(1),L{1}(2),L{1}(3),G{1}(1),G{1}(2),G{1}(3));
    vec1.MaxHeadSize = 0.5;
    vec1.LineWidth = 2.8;

    vec1 = gp.draw('vector rgb:ff0000',G{1}(1),G{1}(2),G{1}(3),X{1}(1),X{1}(2),X{1}(3));
    vec1.MaxHeadSize = 0.5;
    vec1.LineWidth = 2.2;

    end


    if(0)
 %bi2se3

G = {[0 0 0],'$$\Gamma$$'};
Z = {[0 0 0.5]*Kvec,'Z'};
F = {[0.5 0.5 0]*Kvec,'F'};
L = {[0.5 0 0]*Kvec,'L'};

    gmp = gp.draw('point black',G{1}(1),G{1}(2),G{1}(3));
    tg=gp.set_text(gmp,'$$\Gamma$$','Interpreter','Latex');
    tg.FontSize = 16;
    
    %tg.Position = [-0.2 -0.3 -0.4]

    mp = gp.draw('point black',Z{1}(1),Z{1}(2),Z{1}(3));
    tm=gp.set_text(mp,'Z','Interpreter','Latex');
    tm.FontSize = 16;

    kp = gp.draw('point black',F{1}(1),F{1}(2),F{1}(3));
    tk=gp.set_text(kp,'F','Interpreter','Latex');
    tk.FontSize = 16;


    kp = gp.draw('point black',L{1}(1),L{1}(2),L{1}(3));
    tk=gp.set_text(kp,'L','Interpreter','Latex');
    tk.FontSize = 16;


    vec1 = gp.draw('vector rgb:ff0000',G{1}(1),G{1}(2),G{1}(3),Z{1}(1),Z{1}(2),Z{1}(3));
    vec1.MaxHeadSize = 0.5;
    vec1.LineWidth = 2.8;

    vec1 = gp.draw('vector rgb:ff0000',Z{1}(1),Z{1}(2),Z{1}(3),F{1}(1),F{1}(2),F{1}(3));
    vec1.MaxHeadSize = 0.5;
    vec1.LineWidth = 2.2;

    vec1 = gp.draw('vector rgb:ff0000',F{1}(1),F{1}(2),F{1}(3),G{1}(1),G{1}(2),G{1}(3));
    vec1.MaxHeadSize = 0.5;
    vec1.LineWidth = 2.8;

    vec1 = gp.draw('vector rgb:ff0000',G{1}(1),G{1}(2),G{1}(3),L{1}(1),L{1}(2),L{1}(3));
    vec1.MaxHeadSize = 0.5;
    vec1.LineWidth = 2.2;       
    end