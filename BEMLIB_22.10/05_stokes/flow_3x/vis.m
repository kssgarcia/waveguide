close all
clear all

%---
file1 = fopen('flow_3x.net1');
%---
NSG  = fscanf(file1,'%d',[1,1]);

for iseg=1:NSG
  NE(iseg)  = fscanf(file1,'%d',[1,1]);
  for j=1:NE(iseg)
   idle  = fscanf(file1,'%d',[1,1]);
   X(iseg,j) = fscanf(file2,'%f',[1,1]);
   Y(iseg,j) = fscanf(file2,'%f',[1,1]);
  end
end

%---
% plot the geometry
%---

figure(1)
hold on

for iseg=1:NSG
  plot(X(iseg,:),Y(iseg,:),'k')
end



