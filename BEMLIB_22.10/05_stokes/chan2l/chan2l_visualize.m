clear all
close all

%---
file2 = fopen('chan2l.xy');
%---

%---
% read the walls
%---

Ntwo   = fscanf(file2,'%f',[1,1]);
groundhog1   = fscanf(file2,'%f',[4]);
groundhog2   = fscanf(file2,'%f',[4]);

Ntwo   = fscanf(file2,'%f',[1,1]);
groundhog3   = fscanf(file2,'%f',[4]);
groundhog4   = fscanf(file2,'%f',[4]);

%---------
figure(1)
hold on
%---------

%---
% print the walls
%---

plot([groundhog1(2),groundhog2(2)],[groundhog1(3),groundhog2(3)],'r')
plot([groundhog3(2),groundhog4(2)],[groundhog3(3),groundhog4(3)],'r')

axis([-4 4 -2 2])
axis equal
box on
%---------

%---
% loop over profiles
%---

for iloop=1:100

npts  = fscanf(file2,'%f',[1,1]);
if(npts==0) break; end;
whatever   = fscanf(file2,'%f',[1,1]);

nodes  = fscanf(file2,'%f',[6,npts]);

plot(nodes(2,:),nodes(3,:),'k')


end

fclose(file2);
