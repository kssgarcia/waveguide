close all
clear all

%---
% plot a velocity profile from data
% read from file flow_1d.prof
%---

file1 = fopen('flow_1d.prof');

%---
% start reading
%---

Ncl = fscanf(file1,'%f',[1,1]);

for i=1:Ncl
  X0(i) = fscanf(file1,'%f',[1,1]);
  Y0(i) = fscanf(file1,'%f',[1,1]);
  Z0(i) = 0.001;
end

 X0(Ncl+1)=X0(1);
 Y0(Ncl+1)=Y0(1);
 Z0(Ncl+1)=0.001;

 Ndivx  = fscanf(file1,'%f',[1,1]);
 Dx     = fscanf(file1,'%f',[1,1]);
 szx    = fscanf(file1,'%f',[1,1]);
 Ndivy  = fscanf(file1,'%f',[1,1]);
 Dy     = fscanf(file1,'%f',[1,1]);
 szy    = fscanf(file1,'%f',[1,1]);

 for i=1:Ndivx+1
  for j=1:Ndivy+1
    z(i,j) = fscanf(file1,'%f',[1,1]);
  end
 end

%---
% done reading
%---

fclose(file1);

for i=1:Ndivx+1
 x(i) = (i-1)*Dx-szx;
end

for i=1:Ndivy+1
 y(i) = (i-1)*Dy-szy;
end

%---
% plot
%---

figure(1)
hold on

mesh(x,y,z)   % plotting, labelling, and formatting
surf(x,y,z)   % plotting, labelling, and formatting
%waterfall(z)   % plotting, labelling, and formatting

plot3(X0,Y0,Z0,'y')

xlabel('x')
ylabel('y')
zlabel('u')
