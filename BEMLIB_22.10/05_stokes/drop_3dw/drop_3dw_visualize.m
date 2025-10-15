clear all
close all

%---
file2 = fopen('drop_3dw.net');
%---

%---
Npnt   = fscanf(file2,'%f',[1,1]);
Nvert  = fscanf(file2,'%f',[1,1]);
Nface  = fscanf(file2,'%f',[1,1]);

figure(1)
hold on

%----
for iloop=1:100  % over time steps
%----

onemore   = fscanf(file2,'%f',[1,1]);
if(onemore==0) break; end;

clf
%---
vert   = fscanf(file2,'%f',[4,Nvert]);
%---

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
if(colrange<0.001) colrange =1; end 
col = (col-colmin)/colrange;
%col = (col+0.1)/1.5;

coll(:,:,1)=col(:,:);
coll(:,:,2)=col(:,:);
coll(:,:,3)=col(:,:);
%coll = coll+0.25;

%---draw the patches

patch('faces',fac ...
     ,'vertices',vvert ...
     ,'Cdata',col ...
     ,'EdgeColor','k' ...
     ,'FaceColor','interp' ...
     ,'FaceLighting','phong' ...
     ,'BackFaceLighting','lit');
%    ,'FaceColor','y' ...
%    ,'EdgeColor','interp' ...
%light('Position',[1 3 2]);
%light('Position',[-3 -1 3]);
%material dull
%material shiny
%axis vis3d off

Nxy   = fscanf(file2,'%f',[1,1]);
vertxy   = fscanf(file2,'%f',[3,Nxy]);

hold on
plot3(vertxy(1,:),vertxy(3,:),vertxy(2,:),'ro')

%---
% wall
%---

 xw(1)=-1.5; zw(1)=-1.5; yw(1)=0;
 xw(2)= 1.5; zw(2)=-1.5; yw(2)=0;
 xw(3)= 1.5; zw(3)= 1.5; yw(3)=0;
 xw(4)=-1.5; zw(4)= 1.5; yw(4)=0;
 xw(5)=-1.5; zw(5)=-1.5; yw(5)=0;

 patch(xw,zw,yw,yw,'facecolor','y');
 axis([-1.5 1.5 -1.5 1.5 0.0 2.0 ])
 axis equal
 view(30,28)
 box on

pause(1)

end % of iloop

xlabel('x')
ylabel('z')
zlabel('y')

%---
fclose(file2);
