function [s,x,q,beta,alpha,xc,xcl,scl ...
         ,al_scan,xc_scan] ...
 ...
   = static(a,capl,rho1,rho2,kappa,alphain,smax,ndiv)

%-----------
% parameters
%-----------

eps=0.000001;
eps=0.0000001;

tol=0.0000001;

%---
% will scan contact angle space
% with step dal
%---

dal=-0.005;

if(capl<1.5)
 dal=-0.002;
end
if(capl<0.9)
 dal=-0.001;
end
if(capl<0.8)
 dal=-0.0005;
end

Nmax = 2*abs(floor(pi/dal))

%----
% prepare
%---

Drho =rho2-rho1;
Brho =0.5*(rho1+rho2);
rhop = Brho-0.5*kappa*Drho;

capls=capl*capl;  

%----
% flat interface
%---

C(1)=1;
C(2)=0;
C(3)=-3;
C(4)=-4*(Brho-rhop)/Drho;
cosbeta=roots(C);
if(abs(cosbeta(1))<1)
   beta = acos(cosbeta(1));
elseif(abs(cosbeta(2))<1)
   beta = acos(cosbeta(2));
elseif(abs(cosbeta(3))<1)
   beta = acos(cosbeta(3));
end

alpha=beta;
xc=-cos(beta);

%xc=0.05;
%alpha=0.501*pi;
%beta=0.5002*pi;   % initial guess

Ido=1;
Iflag=0;
Icount=0;

%---
while(Ido==1) % over contact angles
%---

Icount=Icount+1;

for iter=1:50

%---
% solve for the floating angle beta
% using Newton's method
%---

  cs  = cos(beta);
  sn  = sin(beta);
  amb = alpha-beta;
  xc  = (2/3*cs*cs*cs +4*(Brho-rhop)/(3*Drho) ...
         +2*capls/a^2*sn*sin(amb))/(1-cs*cs);
  xcl = xc+a*cs;
  scl =    a*sn;
  slope = tan(amb);

  [x,s,q] = flsphere_ode (ndiv,smax,capls,xcl,scl,slope);

  error = x(ndiv+1);

  if(abs(error)<tol) break; end

  beta=beta+eps;

  cs = cos(beta);
  sn = sin(beta);
  amb = alpha-beta;
  xc = (2/3*cs*cs*cs +4*(Brho-rhop)/(3*Drho) ...
       +2*capls/a^2*sn*sin(amb))/(sn*sn);
  xcl = xc+a*cs;
  scl =    a*sn;
  slope = tan(amb);

  [x1,s1,q1] = flsphere_ode(ndiv,smax,capls,xcl,scl,slope);

  error1=x1(ndiv+1);

  beta = beta-eps; % reset

  der = (error1-error)/eps;
  beta = beta-error/der;

end

if(abs(error)>tol) 
   disp "iterations_did_not_converge"
   return
end

%---
% interface profile
%---

  figure(2)
  plot(s,x,'r');
  hold on
%  plot(scl,xcl,'go');
  plot(-s,x,'r');
%  plot(-scl,xcl,'go');

  axis equal
  axis([-2 2 -2 2])
  set(gca,'fontsize',15)
  xlabel('y/a','fontsize',15)
  ylabel('x/a','fontsize',15)

%---
% particle contour
%---

 ncrc = 64;

 for i=1:ncrc+1
  tht = (i-1)*pi*2/ncrc;
  xcrc(i) = xc+a*cos(tht);
  scrc(i) =    a*sin(tht);
 end

% patch(scrc,xcrc,xcrc);
 plot(scrc,xcrc);

 hold off

 pause(0.001)

 if(Iflag==1)
%   pause
   break;
 end

 % alpha_over_pi=alpha/pi
 al_scan(Icount) = alpha/pi;
 xc_scan(Icount) = xc/a;

 alpha=alpha+dal;

if(alpha<0.01*pi)
 dal=abs(dal);
 alpha = alpha+dal;
elseif(alpha>0.99*pi)
 dal = -abs(dal);
 alpha = alpha+dal;
elseif(abs(alpha-alphain)<1.0*abs(dal))  % lock on alphain
 alpha = alphain+0.00;
 Iflag = 1;
end

 if(Icount>Nmax) break; end

%---
end
%---

%---
% done
%---

return
