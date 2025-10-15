%---
file2 = fopen('prtcl_3d.net')
Npnt   = fscanf(file2,'%f',[1,1])
Nvert  = fscanf(file2,'%f',[1,1])
Nface  = fscanf(file2,'%f',[1,1])
vert   = fscanf(file2,'%f',[3,Nvert]);
wall   = fscanf(file2,'%f',[1,1])
fclose(file2)
%---

for i=1:Nvert
save = vert(2,i);
vert(2,i) = vert(3,i);
vert(3,i) = save;
end

Ic=0;
for i=1:Nface
for j=1:Npnt
  Ic=Ic+1;
  fac(j,i) = Ic;
end
end

patch('faces',fac','vertices',vert',...
      'FaceColor','y',...
      'FaceLighting','phong',...
      'BackFaceLighting','lit')
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

xw(1)=-1.5; zw(1)=-1.5; yw(1)=wall;
xw(2)= 1.5; zw(2)=-1.5; yw(2)=wall;
xw(3)= 1.5; zw(3)= 1.5; yw(3)=wall;
xw(4)=-1.5; zw(4)= 1.5; yw(4)=wall;
xw(5)=-1.5; zw(5)=-1.5; yw(5)=wall;
patch(xw,zw,yw,yw);
hold off
