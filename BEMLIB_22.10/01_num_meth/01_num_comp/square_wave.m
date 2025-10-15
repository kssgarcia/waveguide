clear all
close all

%========================================
% FDLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved
%
% This program is to be used only under the
% stipulations of the licensing agreement
%========================================

%===================
% plot a square wave
%===================

figure(1)
hold on
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
set(gca,'fontsize',14)
box on

%---
% plotting points
%---

Nplot=128;

dx = 2.0*pi/Nplot;

for i=1:Nplot+1
 x(i) = (i-1.0)*dx;
end

%---
% nested loops over repeat and Nplot
%---

%---
for repeat=1:6

 if(repeat==1)
  N=1;
 elseif(repeat==2)
  N=5;
 elseif(repeat==3)
  N=9;
 elseif(repeat==4)
  N=17;
 elseif(repeat==5)
  N=33;
 elseif(repeat==6)
  N=65;
 end

 for j=1:Nplot+1
  f(j) = 0.0;
  for i=1:2:N
   f(j) = f(j) + sin(i*x(j))/i;
  end
  f(j) = 4*f(j)/pi;
 end

 if(repeat==1)
  plot(x,f)
 elseif(repeat==2)
  plot(x,f,'r')
 elseif(repeat==3)
  plot(x,f,'y')
 elseif(repeat==4)
  plot(x,f,'c')
 elseif(repeat==5)
  plot(x,f,'k')
 elseif(repeat==6)
  plot(x,f)
 end

end % of repeat
%---
