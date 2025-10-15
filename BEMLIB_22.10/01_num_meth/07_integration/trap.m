close all
clear all

%==========================================
% FDLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement
%==========================================

%======================================
% Integration of a function from a to b
% by the trapezoidal rule
%======================================

%---
% prepare
%---

p = 8;
p = 5;

Nmax = 2^p;

%---
% choose the integral
%---

icase = 13;
icase = 12;
icase = 0;
icase = 11;
icase = 21;
icase = 31;
icase = 72;
icase = 81;

%---
% integration limits
%---

a = 0.0;
b = 1.0;

if(icase==10)
 a = 0.0;
 b = 2*pi;
elseif(icase==11)
 a = 0.0001;
 b = 10.00000;
elseif(icase==12)
 a = 0.0001;
 b = 500.00000;
elseif(icase==13)
 a = 0.0001;
 b = 7.00000;
elseif(icase==21)
 a = 0.0001;
 b = 9.00000;
elseif(icase==31)
 a = 0.0;
 b = 1.30000;
 L = b-a;
elseif(icase==72)
 a = 0.0;
 b = pi;
 alpha = 1.5;
 p72 = 1;
 q72 = 1;
elseif(icase==81)
 a = 0.0;
 b = 1;
end

%---
% base points and function values
%---

h = (b-a)/Nmax;

for i=1:Nmax+1

 x = a+(i-1.0)*h;

 f(i) = cos(pi*x);
 f(i) = exp(-x*x);
 f(i) = exp(-cos(0.5*pi*x)^8);
 f(i) = cos(2*pi*x);
 f(i) = sin(2*pi*x);
 f(i) = x;

 if(icase==10)
  c = 1.4;
  f(i) = (1-c*cos(x))/(1+c^2-2*c*cos(x));
 elseif(icase==11)
  alpha = 1.98;
  f(i) = 4.0*(exp(-x^2)-1.0+x^2)/x^(1+alpha);
 elseif(icase==12)
  alpha = 0.5;
  f(i) = (1-cos(x))/x^(1.0+alpha);
 elseif(icase==13)
  alpha = 1.9;
  f(i) = 2.0*x*exp(-x^2)/x^alpha;
 elseif(icase==21)
  alpha = 1.99;
  q = 1.2;
  f(i) = 2.0*(    exp(-(q-x)^2) ...
             -2.0*exp(-q^2) ...
             +    exp(-(q+x)^2) ...
              +2.0*(1.0-2.0*q^2)*exp(-q^2)*x^2 )/x^(1.0+alpha);
 elseif(icase==31)
  f(i) = exp(-sin(2*pi*x/L)^2);
 elseif(icase==72)
% f(i) = sin(0.5*x)^alpha*sin(p72*x)*sin(q72*x);
  f(i) = 0.5*sin(0.5*x)^alpha*cos(abs(p72-q72)*x);
 elseif(icase==81)
  f(i) = (1-x^2)*cos(0.5*pi*x);
 end

end

%---
% display
%---

figure(1)
hold on;
axis off

%~~~
for j=0:p
%~~~

 N = 2^j;
 fc = Nmax/N;
 h = (b-a)/N;

 pltx(1) = 1;
 plty(1) = 0.1*j;
 if(j==0) plot(pltx,plty,'yo'); end;
 if(j==1) plot(pltx,plty,'bo'); end;
 if(j==2) plot(pltx,plty,'go'); end;
 if(j==3) plot(pltx,plty,'ro'); end;
 if(j==4) plot(pltx,plty,'co'); end;
 if(j==5) plot(pltx,plty,'bo'); end;
 if(j==6) plot(pltx,plty,'ko'); end;

 trapz = 0.5*f(1);

 if(j>0)
  for i=2:N
    k = 1+(i-1)*fc;
    trapz = trapz + f(k);
    pltx(1) = k;
    if(j==1) plot(pltx,plty,'bo'); end;
    if(j==2) plot(pltx,plty,'go'); end;
    if(j==3) plot(pltx,plty,'ro'); end;
    if(j==4) plot(pltx,plty,'co'); end;
    if(j==5) plot(pltx,plty,'bo'); end;
    if(j==6) plot(pltx,plty,'ko'); end;
  end
 end

 trapz = trapz + 0.5*f(Nmax+1);
 trapz = trapz*h;

 if(icase==11)
%   trapz = trapz - 4/(2.0-alpha) * (b^(2.0-alpha)-a^(2.0-alpha));
    trapz = trapz - 4/(2.0-alpha) * (b^(2.0-alpha)-0.0);
    alc = 1.0-alpha;
    trapz = alpha*trapz/(4*gamma(alc)*sin(alc*pi/2));
 end

 if(icase==21)
    trapz = trapz - 4.0*(1.0-2.0*q^2)*exp(-q^2)/(2.0-alpha) * (b^(2.0-alpha)-0.0);
    alc = 1.0-alpha;
    trapz = alpha*trapz/(4*gamma(alc)*sin(alc*pi/2));
 end

 if(icase==31)
    trapz = 2.0*trapz/L;
 end

 Itrap(j+1) = trapz;
 xtrap(j+1) = j;

 pltx(1) = Nmax+1;
 if(j==0) plot(pltx,plty,'yo'); end;
 if(j==1) plot(pltx,plty,'bo'); end;
 if(j==2) plot(pltx,plty,'go'); end;
 if(j==3) plot(pltx,plty,'ro'); end;
 if(j==4) plot(pltx,plty,'co'); end;
 if(j==5) plot(pltx,plty,'bo'); end;
 if(j==6) plot(pltx,plty,'ko'); end;

%~~~
end
%~~~

Itrap'

figure(2)
hold on;
xlabel('log_2N','fontsize',15)
ylabel('log_2(error)','fontsize',15)
set(gca,'fontsize',15)
plot(xtrap(1:p),log2(abs(Itrap(1:p)-Itrap(p+1))),'o:')
box on

figure(3)
plot(f)

if(icase==72)
 n=abs(p72-q72);
 exact = -gamma(-alpha/2+n)*gamma(alpha+1)*sin(0.5*alpha*pi)/gamma(1+alpha/2+n);
 exact = exact/2^(alpha+1);
 exact
end
