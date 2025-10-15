close all
clear all

%---
% plot the wall shear stress from data
% read from file flow_1d.trc
%---

file1 = fopen('flow_1d.trc');

figure(1)
hold on
box on
xlabel('s/a','fontsize',15)
ylabel('\sigma_s','fontsize',15)

%---
% start reading
%---

N = fscanf(file1,'%f',[1,1]);

while (N>0)

for i=1:N
  idle = fscanf(file1,'%f',[1,1]);
  idle = fscanf(file1,'%f',[1,1]);
  idle = fscanf(file1,'%f',[1,1]);
  al(i) = fscanf(file1,'%f',[1,1]);
  sh(i) = fscanf(file1,'%f',[1,1]);
end
plot(al,sh,'k-')

N = fscanf(file1,'%f',[1,1]);

end


