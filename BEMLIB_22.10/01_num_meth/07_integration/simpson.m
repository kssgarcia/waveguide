close all
clear all

%-----------------------------------
% One-dimensional integration from a to b
% by the Simpson rule
%
% will compute the integral for a
% cascace of divisions up to Nmax
%-----------------------------------

p = 5;
p = 6;

Nmax=2^p;

a = 0.0;
b = 1.0;

h = (b-a)/Nmax;

%---
% function values
%---

for i=1:Nmax+1
 x = a+(i-1.0)*h;
 f(i) = cos(pi*x);
 f(i) = cos(2*pi*x);
 f(i) = sin(2*pi*x);
 f(i) = 1.0;
 f(i) = x;
 f(i) = x^2;
 f(i) = exp(-x*x);
 f(i) = exp(-cos(0.5*pi*x)^8);
 f(i) = x^1.5;
end

hold on;

%----
% loop over grids
%----

%~~~
for j=1:p
%~~~

 N=2^j;
 fc = Nmax/N;
 h = (b-a)/N;

  pltx(1) = 1;
  plty(1) = 0.2*j;

  if(j==1) plot(pltx,plty,'bo'); end;
  if(j==2) plot(pltx,plty,'go'); end;
  if(j==3) plot(pltx,plty,'ro'); end;
  if(j==4) plot(pltx,plty,'co'); end;
  if(j==5) plot(pltx,plty,'bo'); end;
  if(j==6) plot(pltx,plty,'ko'); end;

 Ismp=f(1);
 weight = 4.0;

 for i=2:N
   k=1+(i-1)*fc;
   Ismp = Ismp+weight*f(k);
   if(weight>3.99)
     weight =2.0;
   else
     weight = 4.0;
   end
   pltx(1)=k;
   if(j==1) plot(pltx,plty,'bo'); end;
   if(j==2) plot(pltx,plty,'go'); end;
   if(j==3) plot(pltx,plty,'ro'); end;
   if(j==4) plot(pltx,plty,'co'); end;
   if(j==5) plot(pltx,plty,'bo'); end;
   if(j==6) plot(pltx,plty,'ko'); end;
 end

 Ismp = Ismp+f(Nmax+1);
 Ismp = Ismp*h/3.0;

 Isimp(j) = Ismp;
 xsimp(j)=j;

 pltx(1) = Nmax+1;

 if(j==1) plot(pltx,plty,'bo'); end;
 if(j==2) plot(pltx,plty,'go'); end;
 if(j==3) plot(pltx,plty,'ro'); end;
 if(j==4) plot(pltx,plty,'co'); end;
 if(j==5) plot(pltx,plty,'bo'); end;
 if(j==6) plot(pltx,plty,'ko'); end;

%~~~
end
%~~~

Isimp'

axis off

figure
plot(xsimp(1:p-1),log2(abs(Isimp(1:p-1)-Isimp(p))),'o:')
xlabel('log_2N','fontsize',15)
ylabel('log_2(error)','fontsize',15)
set(gca,'fontsize',15)
