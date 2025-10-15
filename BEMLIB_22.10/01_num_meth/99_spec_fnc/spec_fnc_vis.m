close all
clear all
%----

file2 = fopen('spec_fnc.out')

%---
% prepare to plot
%---

figure(1)
hold on
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
box on
axis([0 8 -10 10])

N = 1

ipass = 0;

%---
while(N~=0)
%---

N  = fscanf(file2,'%d',[1,1]);

for i=1:N
 idle  = fscanf(file2,'%d',[1,1]);
 x(i)  = fscanf(file2,'%f',[1,1]);
 y(i)  = fscanf(file2,'%f',[1,1]);
end

ipass = ipass+1;

if(ipass==1)
 plot(x,y,'k')
elseif(ipass==2)
 plot(x,y,'k--')
end

%---
end
%---
