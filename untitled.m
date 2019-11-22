close all;
clear;
clc;
f = @(x,y,z) x.^2 + y.^2 + z.^2 - 5;
interval = [-5 5 -5 5 -2 5];
%fimplicit3(f,interval)

a = 1e-10;
from = -2*pi./a;
to = -from;
len = 100;
k = linspace(from,to,len);
[kx,ky,kz] = meshgrid(k);
Eo = 0;
t = 0.5;
En = -Eo - 4.*t.*(cos(ky.*a/2).*cos(kz.*a/2)+cos(kx.*a/2).*cos(kz.*a/2)+cos(kx.*a/2).*cos(ky.*a/2));
Ef1 = -Eo - 8*t;

E = @(kx,ky,kz) -Eo - 4.*t.*(cos(ky.*a/2).*cos(kz.*a/2)+cos(kx.*a/2).*cos(kz.*a/2)+cos(kx.*a/2).*cos(ky.*a/2)) + Eo - 1*t;
int = [from to from to from to];
fimplicit3(E,int);

%contour3(kx(:,:,1),ky(:,:,1),E(:,:,1))
%surf(kx(:,:,1),ky(:,:,1),E(:,:,1))

[kx,ky] = meshgrid(k);
Esq = -Eo - 2.*t.*(cos(kx.*a)+cos(ky.*a));
figure();
contour(kx,ky,Esq);



       % function set_orbital(self,varargin)
          
       %    orbs = struct('name','s','l',0,'orbital_number',1);
       %    orbp = struct('name','p','l',1,'orbital_number',3);
       %    orbd = struct('name','d','l',2,'orbital_number',5);
       %    orbf = struct('name','f','l',3,'orbital_number',7);          

       %    possible_orbitals = {orbs,orbp,orbd,orbf};

       %    for i = 1:size(varargin,2)
       %      if(isstring(varargin{i}) | ischar(varargin{i}))
       %        o = '';
       %        for orb_iter = 1:size(possible_orbitals,2)
       %          o = possible_orbitals{orb_iter};
       %          user_orb = varargin{i};
       %          n = str2num(user_orb(1));
       %          orb_name = user_orb(2:end);
       %          %fprintf("user %s, us %s ,cmp = %d\n",orb_name,o.name,strcmp(orb_name,o.name));
       %          if(strcmp(orb_name,o.name))
       %            break;
       %          end
       %        end
              
       %        if(isempty(o) == 0)
       %          self.no_orbital = self.no_orbital + o.orbital_number;
       %          self.orbitals{end+1} = user_orb;
       %        else
       %          fprintf('Entered orbital %s\n',orb_name)
       %          error('Undefined orbital it must be s,p,d or f');
       %        end
              
       %      end
       %    end
       %    %self.no_orbital = size(self.orbitals,2);
       % end