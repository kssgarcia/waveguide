clear all
close all

%---
file2 = fopen('caps_2d.xy');
%---

%---------
figure(1)
hold on
axis([-1.5 1.5 -1.5 1.5])
axis equal
box on
%---------

%---
% loop over profiles
%---

for iloop=1:200

npts  = fscanf(file2,'%f',[1,1]);
if(npts==0) break; end;
whatever   = fscanf(file2,'%f',[1,1]);
whatever   = fscanf(file2,'%f',[1,1]);

nodes  = fscanf(file2,'%f',[6,npts]);

plot(nodes(2,:),nodes(3,:),'.k-')

end

fclose(file2);
