close all
clear all
hold on;

%=================================
% solution of system of two ODEs
% using the Pozrikidis/Prony method
%
% the starting series is generated
% analytically
%=================================

%---
Dt = 0.1;
Dt = 0.2;
%---

%---
% generate the starting points
%---

Iselect=1;

if(Iselect==1) %--- linear
  b=0.3;
  x0 = 1.0;
  x1 = b+(x0-b)*exp(Dt);
  x2 = b+(x0-b)*exp(2*Dt);
  Nmax=100;
elseif(Iselect==2) %--- cosine
  x0 = 1.0;
  x1 = cos(Dt);
  x2 = cos(2*Dt);
  Nmax=10;
elseif(Iselect==3) %--- sine
  x0 = 0.0;
  x1 = sin(Dt);
  x2 = sin(2*Dt);
  Nmax=200;
end

plot(0,x0,'o')
plot(Dt,x1,'o')
plot(2*Dt,x2,'o')

a2 = -(x2-x1)/(x1-x0);

X = (x1+a2*x0)/(1+a2);

%---
% now continue
%---

%---
x(2)=x2;
%---

for j=3:Nmax

 x(j)=(1+a2)*X-a2*x(j-1);

 time = j*Dt;
 if(Iselect==1)
  exact = b+(x0-b)*exp(time);
 elseif(Iselect==2)
  exact = cos(time);
 elseif(Iselect==3)
  exact = sin(time);
 end
 plot(time,x(j),'o')
 plot(time,exact,'r+')
 if(time>=10) break; end

end
%---
