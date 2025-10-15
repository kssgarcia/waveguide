clear all
close all

%---
file2 = fopen('bump_3d.net');
%---

%---
Npnt   = fscanf(file2,'%f',[1,1]);
Nvert  = fscanf(file2,'%f',[1,1]);
Nface  = fscanf(file2,'%f',[1,1]);
vert   = fscanf(file2,'%f',[4,Nvert]);
wall   = fscanf(file2,'%f',[1,1]);
%---

%---vertices

for i=1:Nvert
  vvert(i,1) = vert(1,i);
  vvert(i,2) = vert(3,i);
  vvert(i,3) = vert(2,i);
end

%---color

Ic=0;
for i=1:Nface
 for j=1:Npnt
  Ic=Ic+1;
  col(j,i) = vert(4,Ic);
 end
end

%---faces
% vertex connection defining each face (connectivity matrix)

Ic=0;
for i=1:Nface
for j=1:Npnt
  Ic=Ic+1;
  fac(i,j) = Ic;
end
end

%---scale the color

colmax = -200.0;
colmin =  200.0;

for i=1:Nface
for j=1:Npnt
   if(col(j,i)>colmax) colmax = col(j,i); end
   if(col(j,i)<colmin) colmin = col(j,i); end
end
end
colrange = colmax-colmin
col=(col-colmin)/colrange;

%---draw the patches

patch('faces',fac,...
     'vertices',vvert,...
       'Cdata',col,...
      'EdgeColor','r',...
     'FaceColor','interp',...
      'FaceLighting','phong',...
      'BackFaceLighting','lit')
%'FaceColor','y',...
%light('Position',[1 3 2]);
%light('Position',[-3 -1 3]);
%material dull
%material shiny
%axis vis3d off
axis([-1.5 1.5 -1.5 1.5 -0.5 2.5 ])
%view(45,34)
xlabel('x')
ylabel('z')
zlabel('y')

%---
% wall
%---

xw(1)=-1.5; zw(1)=-1.5; yw(1)=wall;
xw(2)= 1.5; zw(2)=-1.5; yw(2)=wall;
xw(3)= 1.5; zw(3)= 1.5; yw(3)=wall;
xw(4)=-1.5; zw(4)= 1.5; yw(4)=wall;
xw(5)=-1.5; zw(5)=-1.5; yw(5)=wall;
patch(xw,zw,yw,yw,'facecolor','y');
hold on

%---
% collocation points
%---

ncl  = fscanf(file2,'%f',[1,1]);
points = fscanf(file2,'%f',[3,ncl]);

for i=1:ncl
 xp(i) = 1.02*points(1,i);
 yp(i) = 1.02*points(2,i);
 zp(i) = 1.02*points(3,i);
end

%plot3(xp,zp,yp,'o')

break

%---
% streamlines
%---

nev = fscanf(file2,'%f',[1,1]);

while (nev>0)

streamline  = fscanf(file2,'%f',[3,nev]);

for i=1:nev
 xx(i) = streamline(1,i);
 yy(i) = streamline(2,i);
 zz(i) = streamline(3,i);
end
plot3(xx,zz,yy);
clear xx yy zz
nev = fscanf(file2,'%f',[1,1]);

end

fclose(file2);
