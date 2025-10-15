close all
clear all

%=======================
% visualize a net of 6-node triangles
%
% data are read from file: trgl6.net
%=======================

file2 = fopen('trgl6.net');

Npnt   = fscanf(file2,'%f',[1,1])
Nvert  = fscanf(file2,'%f',[1,1])
Nface  = fscanf(file2,'%f',[1,1])
vert   = fscanf(file2,'%f',[4,Nvert]);

%---
%Ic = 0
%Jc = 0
%Np = fscanf(file2,'%f',[1,1])
%Npnt = Np
%for iter=1,5000
%Jc = Jc+1
% for i=1,Np
%   Ic=Ic+1
%   vert(1:3,Ic) = fscanf(file2,'%f',[3,1])
%end 
%Np = fscanf(file2,'%f',[1,1]);
% If Np < 1
%  break
% end
%end
%Nvert = Ic
%Nface = Jc
%---

fclose(file2)

%---vertices

for i=1:Nvert
  vvert(i,1) = vert(1,i);
  vvert(i,2) = vert(3,i);
  vvert(i,3) = vert(2,i);
end

%---color

Ic = 0;

for i=1:Nface
 for j=1:Npnt
   Ic=Ic+1;
   col(j,i) = vert(4,Ic);
 end
end

%---faces

% vertex connection defining each face (connectivity matrix)

Ic = 0;

for i=1:Nface
 for j=1:Npnt
  Ic = Ic+1;
  fac(i,j) = Ic;
 end
end

patch('faces',fac ...
     ,'vertices',vvert ...
     ,'Cdata',col ...
     ,'FaceColor','y' ...
     ,'FaceLighting','phong' ...
     ,'BackFaceLighting','lit')
%light('Position',[1 3 2]);
%light('Position',[-3 -1 3]);
%material dull
%material shiny
%axis vis3d off
%view(45,34)
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
zlabel('z','fontsize',15)
set(gca,'fontsize',15)
box on
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5 ])
axis equal

hold off
