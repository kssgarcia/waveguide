%close all
clear all

%---
file2 = fopen('trgl6.net');
%---

%---
Npnt   = fscanf(file2,'%d',[1,1]);
Nvert  = fscanf(file2,'%d',[1,1]);
Nface  = fscanf(file2,'%d',[1,1]);
% strain  = fscanf(file2,'%f',[1,1]);
vert   = fscanf(file2,'%f',[4,Nvert]);
%---

a11 = 1; a12 = 0;
a21 = strain; a22 = 1;


Ic=0;
for i=1:Nface
hold on
 for j=1:Npnt
  Ic=Ic+1;
  xplot(j)=vert(1,Ic);
  yplot(j)=vert(2,Ic);
  zplot(j)=vert(3,Ic);
 end
 plot3(xplot,yplot,zplot,'-')
%plot3(xplot(1),yplot(1),zplot(1),'o')
%plot3(xplot(2),yplot(2),zplot(2),'+')
end

Np = 1;

for ip=-Np:Np
for jp=-Np:Np

%---vertices

for i=1:Nvert
  vvert(i,1) = vert(1,i)+ip*a11+jp*a21;
  vvert(i,2) = vert(2,i)+ip*a12+jp*a22;
  vvert(i,3) = vert(3,i);
end

%---color

Ic=0;
for i=1:Nface
 for j=1:Npnt
  Ic=Ic+1;
  col(j,i) = vert(4,Ic);
 end
end

%---
% faces
%---
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

colrange = colmax-colmin;
col = (col-colmin)/colrange;

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

end
end

%axis([-1.5 1.5 -1.5 1.5 -0.5 2.5 ])
%view(45,34)
xlabel('x')
ylabel('y')
zlabel('z')
axis square
axis equal

fclose(file2);

%---

figure(2)
hold on
Ic=0;
for i=1:Nface
hold on
 for j=1:Npnt
  Ic=Ic+1;
  xplot(j)=vert(1,Ic);
  yplot(j)=vert(2,Ic);
  zplot(j)=vert(4,Ic);
 end
 plot3(xplot,yplot,zplot,'-')
end

