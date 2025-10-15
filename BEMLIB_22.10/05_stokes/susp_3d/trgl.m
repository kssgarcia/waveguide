file2 = fopen('trgl.net')
Npnt  = fscanf(file2,'%f',[1,1])
Nvert = fscanf(file2,'%f',[1,1])
Nface = fscanf(file2,'%f',[1,1])
vert  = fscanf(file2,'%f',[3,Nvert]);

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

Ic=0;
for i=1:Nface
for j=1:Npnt
  Ic=Ic+1;
  fac(j,i) = Ic;
end
end

%patch('faces',fac','vertices',vert',...
%      'FaceColor','y',...
%      'FaceLighting','phong',...
%      'BackFaceLighting','lit')
patch('faces',fac','vertices',vert',...
      'FaceColor','w',...
      'FaceLighting','phong',...
      'BackFaceLighting','lit')
%light('Position',[1 3 2]);
%light('Position',[-3 -1 3]);
%material dull
%material shiny
%axis vis3d off
axis([-1.0 1.0 -1.0 1.0 -1.0 1.0 ])
%view(45,34)
xlabel('x')
ylabel('y')
zlabel('z')

hold on

break

Nsp = fscanf(file2,'%f',[1,1])
col  = fscanf(file2,'%f',[3,Nsp]);

for i=1:Nsp
 plot3(col(1,i),col(2,i),col(3,i),'b+')
  hold on
end

fclose(file2)

