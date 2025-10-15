close all
clear all

%----------------------------------------------
% Driver for flow
% due to a periodic array of point-force rings
% inside a cylinder of radius sc
%----------------------------------------------

figure(1)
hold on

sc = 1.0;  % tube radius
RL = 1.0;  % period
x0 = 0.0;  % location of the ring
s0 = 0.5;  % radius of the ring

Np = 2;

Nprof=16;
dprof = (sc-0.01)/Nprof;

%--------
for irepeat1=1:2
%--------

if(irepeat1==1)
 Nsum=2;
elseif(irepeat1==2)
 Nsum=4;
end

%--------
for irepeat=1:2
%--------

if(irepeat==1)
 x = 0.1;
elseif(irepeat==2)
 x = 0.5;
end

%--------
for j=1:Nprof+1
%--------

 s = 0.01+(j-1.0)*dprof;

  [gxx,gxs  ...
  ,gsx,gss] ...
...
      = sgf_ax_1p_ct  ...
...
  (x,s     ...
  ,x0,s0   ...
  ,sc      ...
  ,RL      ...
  ,Nsum,Np ...
  );

  radial(j) = s;
  velx(j) = gxx;

end
%------------

if(irepeat1==1)
 plot(velx,radial,'-o')
 plot(velx,-radial,'-o')
elseif(irepeat1==2)
 plot(velx,radial,'-ro')
 plot(velx,-radial,'-ro')
end

%---
end
%---

%---
end
%---

figure(1)
plot([-2.0,2.0],[sc,sc],'k-')
plot([-2.0,2.0],[-sc,-sc],'k-')
xlabel('u','fontsize',15)
ylabel('\sigma','fontsize',15)
axis([-2.0 2.0, -1.1*sc,1.1*sc])
set(gca,'fontsize',15)

