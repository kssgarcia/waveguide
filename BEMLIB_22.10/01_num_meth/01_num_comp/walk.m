clear all
close all

%=====================================
% Simulation of discrete random walks
%
% At the initial instant, the walkers
% are stacked in two columns at x=+-1/2
%
% Left and right walkers are treated separately
%
% M:       number of walkers
% notches: number of notches
% Nstep:   number of steps
% q:       left probability
%=====================================

%---
% settings
%---

M = 2*512;
notches = 512;

Nstep = 50;
Nstep = 10;
Nstep = 200;

q = 0.75;
q = 1.00;
q = 0.5;

%---
% prepare
%---

M2 = 2*M;

%-----------------------
% right and left notches
% and initial populations
%-----------------------

for i=1:notches
  xr(i) = i-0.5;
  kr(i) = 0;
  xl(i) =-i+0.5;
  kl(i) = 0;
end

for i=1:2*notches
  xgauss(i) = -notches+i-0.5;
  ygauss(i) =  0.0;
end

%-----------
% initialize
%-----------

kr(1) = M;
kl(1) = M;

for i=1:M2
  yplot_r(i)=0;
  xplot_r(i)=0;
  yplot_l(i)=0;
  xplot_l(i)=0;
end

%===============
for step=1:Nstep
%===============

%---
% plotting
%---

figure(1)

%clear xplot_r
%clear yplot_r
%clear xplot_l
%clear yplot_l

Icr = 0; % counter

for i=1:notches

 if(kr(i)>0)
  for j=1:kr(i)
   Icr = Icr+1;
   xplot_r(Icr) = xr(i);
   yplot_r(Icr) = 1.0*j;   % stack particles vertically
  end
 end

end

Icl = 0;

for i=1:notches

 if(kl(i)>0)
  for j=1:kl(i)
   Icl = Icl+1;
   xplot_l(Icl) = xl(i);
   yplot_l(Icl) = 1.0*j;   % stack particles vertically
  end
 end

end

%---
% animation
%---

if(step==1)
  Handle1 = plot(xplot_r,yplot_r/M2,'k.');
  hold on
  Handle2 = plot(xplot_l,yplot_l/M2,'k.');
  Handle3 = plot(xgauss,ygauss,'r');
  set(Handle2,'EraseMode','xor')
  set(Handle1,'EraseMode','xor')
  set(Handle3,'EraseMode','xor')
  xlabel('i','fontsize',15)
  ylabel('k/(2M)','fontsize',15)
  set(gca,'fontsize',15)
end

%=========
% stepping
%=========

krnew = kr;
klnew = kl;

avr = 0.0;  % average
var = 0.0;  % variance

if(kr(1)>0)

  for j=1:kr(1)
   krnew(1) = krnew(1)-1;
   if(rand>q)
    klnew(1) = klnew(1)+1;
   else
    krnew(2) = krnew(2)+1;
   end
  end

  avr = avr+kr(1)*xr(1);
  var = var+kr(1)*xr(1)^2;

end

if(kl(1)>0)

  for j=1:kl(1)
   klnew(1) = klnew(1)-1;
   if(rand>q)
    klnew(2) = klnew(2)+1;
   else
    krnew(1) = krnew(1)+1;
   end
  end

  avr = avr+kl(1)*xl(1);
  var = var+kl(1)*xl(1)^2;

end

%---
for i=2:notches

 if(kr(i)>0)

   for j=1:kr(i)
    krnew(i) = krnew(i)-1;
    if(rand>q)
     krnew(i-1) = krnew(i-1)+1;
    else
     krnew(i+1) = krnew(i+1)+1;
    end
  end

  avr = avr+kr(i)*xr(i);
  var = var+kr(i)*xr(i)^2;

 end

 if(kl(i)>0)

  for j=1:kl(i)

   klnew(i) = klnew(i)-1;
   if(rand>q)
    klnew(i+1) = klnew(i+1)+1;
   else
    klnew(i-1) = klnew(i-1)+1;
   end
  end

  avr = avr+kl(i)*xl(i);
  var = var+kl(i)*xl(i)^2;

 end

 end  % over notches
%---

time(step) = step-1;
mean(step) = avr/M2;
variance(step) = -mean(step)^2+var/M2;

s = sqrt(variance(step));

for i=1:2*notches
  ygauss(i) = 1/sqrt(2*pi)/s*exp(-0.5*xgauss(i)^2/s^2);
end


kr = krnew;
kl = klnew;

%total=0;
%for i=1:notches
% total = total+kr(i)+kl(i);
%end
%total

%axis([-step, step, 0, max(max(kr),max(kl))])
% axis([-step, step, 0, 1])
% axis([-step, step, 0, max(max(kr),max(kl))/M2])
 axis([-100, 100, 0,  max(max(kr),max(kl))/M2])
 set(Handle1,'XData',xplot_r,'YData',yplot_r/M2)
 set(Handle2,'XData',xplot_l,'YData',yplot_l/M2);
 set(Handle3,'XData',xgauss,'YData',ygauss);
 box on
 drawnow
 pause(0.5)

 for i=1:M2
  xplot_r(i) = 0;
  yplot_r(i) = 0;
  xplot_l(i) = 0;
  yplot_l(i) = 0;
 end


%---
end
%---

figure(2)
plot(time,mean,'k')
xlabel('N','fontsize',15)
ylabel('Mean','fontsize',15)
set(gca,'fontsize',15)

figure(3)
plot(time,variance,'k')
xlabel('N','fontsize',15)
ylabel('Variance','fontsize',15)
set(gca,'fontsize',15)
