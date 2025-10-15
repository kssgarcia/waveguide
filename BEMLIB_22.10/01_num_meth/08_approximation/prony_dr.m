clear all

%=====================================
% Prony fitting of a signal consisting
% of a sum of complex exponentials
%=====================================

%---
% data
%---

Np=32;

%---
% coefficients
%---

c0=0.2;
K=4;
j=1;
sig(j,1)=0.35; sig(j,2)=-0.34; ac(j)=4.56; as(j)=3.491;
j=2;
sig(j,1)=1.45; sig(j,2)=-0.45; ac(j)=0.56; as(j)=2.491;
j=3;
sig(j,1)=1.45; sig(j,2)=-0.45; ac(j)=0.56; as(j)=2.491;
sig(j,1)=1.15; sig(j,2)= 0.90; ac(j)=2.56; as(j)=1.491;
j=4;
sig(j,1)=-0.45; sig(j,2)=1.45; ac(j)=0.06; as(j)=1.491;


Dt = 0.1/min(abs(sig(:,2)));
Dt = 0.2;

%---
% generate the signal
%---

for i=1:Np
   t = (i-1.0)*Dt;
   f(i) = c0;
   for j=1:K
     arg = sig(j,1)*t;
     comp = (ac(j)*cos(arg)+as(j)*sin(arg)) * exp(sig(j,2)*t);
%      if(j==3)
%       comp=comp*t;
%      end
     f(i) = f(i) + comp;
    end
end

[sigma] = prony (Np,f,Dt,K);

sig
sigma
