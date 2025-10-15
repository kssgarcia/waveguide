clear all
close all

%----------------
% display the triangulation of a section of a sphere (hsph)
%
% data are read from tile: trgl6_hsph.net
%----------------

%---
file2 = fopen('trgl6_hsph.net')
%---

Npnt   = fscanf(file2,'%f',[1,1])
Nvert  = fscanf(file2,'%f',[1,1])
Nface  = fscanf(file2,'%f',[1,1])

vert   = fscanf(file2,'%f',[4,Nvert]);

fclose(file2)
%---

%---vertices

for i=1:Nvert
  vvert(i,1) = vert(3,i);
  vvert(i,2) = vert(1,i);
  vvert(i,3) = vert(2,i);
end

%---color

Ic = 0;

for i=1:Nface
 for j=1:Npnt
  Ic = Ic+1;
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


figure(1)
hold on

patch('faces',fac ...
     ,'vertices',vvert ...
     ,'Cdata',col ...
     ,'FaceColor','y' ...
     ,'FaceLighting','phong',...
      'BackFaceLighting','lit')
%light('Position',[1 3 2]);
%light('Position',[-3 -1 3]);
%material dull
%material shiny
%axis vis3d off
view(125,22)

xw(1)=-1.5;yw(1)=-1.5;zw(1)=0.0;
xw(2)= 1.5;yw(2)=-1.5;zw(2)=0.0;
xw(3)= 1.5;yw(3)= 1.5;zw(3)=0.0;
xw(4)=-1.5;yw(4)= 1.5;zw(4)=0.0;

patch(xw,yw,zw ...
     ,'FaceColor','r')

xlabel('z','fontsize',15)
ylabel('x','fontsize',15)
zlabel('y','fontsize',15)
set(gca,'fontsize',15)
axis([-1.5 1.5 -1.5 1.5 -1.5 2.5 ])
axis equal
box on


