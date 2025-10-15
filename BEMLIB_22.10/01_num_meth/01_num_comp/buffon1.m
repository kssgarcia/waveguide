clear all
close all

%=========================
% buffon needle simulation
%=========================

figure(1)
hold on
xlabel('log(N)','fontsize',15)
ylabel('log(|\pi_{MC}-\pi|)','fontsize',15)
set(gca,'fontsize',15)
box

%---
% line of slope -1/2
%---

xplot(1)= 0.0; yplot(1)= 2.0;
xplot(2)=10.0; yplot(2)=-3.0;

plot(xplot,yplot,'r');

%---
% settings
%---

a=1.0;
b=0.9;
boa=b/a;

%---
% Monte-Carlo
%---

Nmax=1024*2*2*2*2*2;

M=0;
jplot=32;

for N=1:Nmax

  r1 = rand;
  r2 = rand;
  theta = pi*r2;
  sn = sin(theta);

  if(r1<boa*sn)
    M=M+1;
  end

  P=M/N;

  if(P>0)
   pi_mc = 2.0*boa/P;
  else
   pi_mc=0.0;
  end

  if(jplot==32)
  plot(log(N),log(abs(pi_mc-pi)),'.')
  jplot=0;
  end

  jplot=jplot+1;

end
