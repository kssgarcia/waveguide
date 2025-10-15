%----
close all
clear all
%----

%============
% visualization
%============

file2 = fopen('prd_3d_pr.net')

%---

N   = fscanf(file2,'%f',[1,1])
nodes  = fscanf(file2,'%f',[4,N]);

figure(1)
hold on
axis equal
plot3(nodes(2,:),nodes(3,:),nodes(4,:),'ko-')
box on
view(47,28)

Nnew   = fscanf(file2,'%f',[1,1])
nodesnew  = fscanf(file2,'%f',[4,Nnew]);

plot3(nodesnew(2,:),nodesnew(3,:),nodesnew(4,:),'r.')
%---

fclose(file2)

