%---
clear all
close all
%---

%===========================
% motion of a viscous drop
%===========================

%---

  Iflow = 3;   % inside a cylinder
  Iflow = 1;   % free space
  Iflow = 2;   % above a wall
%---

%---
  mu = 1.0;  % viscosity
  lambda = 1.0;
  gac =1.0;   % acceleration of gravity
  rhod=1.0;  % density of the drop
  rhoa=0.0;  % density of the ambient fluid
  a = 1.0;   % drop radius

  gamma=2.00; % surface tension
  gamma=1.00; % surface tension
  gamma=0.25; % surface tension
  gamma=0.20; % surface tension
  gamma=0.00; % surface tension

  wall = 0.0;

  if(Iflow==1)
   xcnt = 0.0;
   xcnt = 1.1;
  elseif(Iflow==2)
   xcnt = 2.0;
  elseif(Iflow==3)
   xcnt = 0.0;
  end

  sc = 1.2;
  RL = 4.0;
  Nsum = 2;
  Np = 2;

%---
% prepare
%---

  terminal = - (rhod-rhoa)/mu * 2/3 *(lambda+1)/(3*lambda+2) *gac;

  VD = 4.0*pi*a^3/3;   % drop volume

  capls = gamma/(gac*abs(rhod-rhoa));  % square of the capillary number

%---
% interface perturbation amplitude
%---

 amp=-0.2*a;
 amp= 0.0*a;
 amp= 0.2*a;

%--------
% interfacial segments
%-------

 NSG=08;
 NSG=256;
 NSG=128;
 NSG=64;
 NSG=96;
 NSG=32;
 NSG=16;

%--------
% Gauss-Legendre quadrature
%-------

 NGL=2;
 NGL=12;
 NGL=6;

%---
% time stepping
%---

IRK = 2; % RK2
IRK = 1; % Euler

Dt=5.0; 

Dt=2.5;
Dt=1.0;
Dt=2.5;
Dt=0.075;
Dt=0.050;
Dt=0.05;
Dt=0.020;
Dt=0.1;
Dt=0.5;


move=1; % normal velocity
move=0; % total velocity

Nstep=2;
Nstep=40;
Nstep=1;
Nstep=10;
Nstep=30;
Nstep=50;
Nstep=5;
Nstep=20;
Nstep=60;
Nstep=180;
Nstep=3600;
Nstep=2400;
Nstep=1200;
Nstep=600;
Nstep=1800;
Nstep=400;

%---
% integration over the interface
%---

 Ichoose=1; % lines
 Ichoose=2; % splines

%---
% point redistribution
%---

Ich1=0;
Ich1=1;
thmax=0.25*pi;
thmax=0.030*pi;
thmax=0.125*pi;

Ich2=0;
Ich2=1;
spmax=1.0*pi*a/NSG;
spmax=2.0*pi*a/NSG;

Ich3=0;
Ich3=1;
spmin=spmax/5.0;
spmin=spmax/2.0;

%-------------
% initial shape
%-------------

   for i=1:NSG+1
    angle = (i-1.0)*pi/NSG;
    rad = a + amp*cos(2*angle);
    XI(i) = xcnt + rad*cos(angle);
    YI(i) =        rad*sin(angle);
   end

%-------------
% plotting
%-------------

 figure(1)
 hold on
 axis equal
 plot( YI,-XI,'-k.')
 plot(-YI,-XI,'r:')
 box on

%----------------
% open data files
%----------------

file1 = fopen('drop_ax.out','w');

fprintf(file1,'%12.8f ',VD);
fprintf(file1,'\n');
fprintf(file1,'\n');
fprintf(file1,'%12.8f ',amp);
fprintf(file1,'\n');
fprintf(file1,'%12.8f ',thmax);
fprintf(file1,'\n');
fprintf(file1,'%12.8f ',spmax);
fprintf(file1,'\n');
fprintf(file1,'%12.8f ',spmin);
fprintf(file1,'\n');
fprintf(file1,'\n');

%-----------
% initialize
%-----------

time(1)=0.0;

%================
for step=1:Nstep
%================

step

%---
% animation
%---

clear XXI YYI

for i=1:NSG+1
  XXI(i)= XI(i);
  YYI(i)= YI(i);
  if(Iflow==1)
    XXI(i) = XXI(i)+terminal*time(step);
  end
end

for i=1:NSG
  XXI(NSG+1+i) =  XI(NSG+1-i);
  YYI(NSG+1+i) = -YI(NSG+1-i);
  if(Iflow==1)
    XXI(NSG+1+i) = XXI(NSG+1+i)+terminal*time(step);
  end
end

if(step==1)
  figure(2)
  Handle1 = patch( YYI, XXI, 'y');
  hold on
  Handle2 = plot( YYI, XXI, '.-');
  axis equal
  set(gca,'fontsize',12)
  xlabel('x/a','fontsize',12)
  ylabel('y/a','fontsize',12)
%  axis([-1.2  1.2 -1.6 1.6])
   axis([-1.5 1.5 -8  1])
   if(Iflow==2) 
      patch([-2.0 2.0 2.0 -2.0],[0.0 0.0 -0.1 -0.2],'r')
      axis([-2.0 2.0 -0.1 4])
   end
   if(Iflow==3)
      plot( YYI, XXI+RL, '.-')
      plot( YYI, XXI-RL, '.-')
      plot( YYI, XXI-2*RL, '.-')
   end
  box on
else
  set(Handle1,'XData',[YYI] ...
             ,'YData',[XXI]);
  set(Handle2,'XData',[YYI] ...
             ,'YData',[XXI]);
  drawnow
   if(Iflow==3)
      plot( YYI, XXI+RL, '.-')
      plot( YYI, XXI-RL, '.-')
      plot( YYI, XXI-2*RL, '.-')
   end
%  axis([-1.2  1.2 -1.6 1.6])
end

%if(step==1 | step==80 | step==160 | step==240 | step==320)
  figure(5)
  hold on
  axis equal
% patch( YYI, XXI, 'y');
%  plot( YYI, XXI, 'k.-');
  plot( YYI, XXI, 'k-');
  axis([-1.5 1.5 -11  1])
  axis([-1.5 1.5 -9  1.5])
  if(Iflow==3) plot( YYI, XXI+RL, '.-'); end
  if(Iflow==3) plot( YYI, XXI-RL, '.-'); end
  if(Iflow==3) plot( YYI, XXI-2*RL, '.-'); end
  if(Iflow==2) axis([-2.0 2.0 -0.1 4]); end
  if(Iflow==2) patch([-2.0 2.0 2.0 -2.0],[0.0 0.0 -0.2 -0.2],'r'); end
  box on
%end

%---
% point redistribution
%---

Ido=1;
Iloop=0;

%--
while(Ido==1)
%--

  Iloop=Iloop+1;

  [NSG,XI,YI ...
  ,vnx,vny,crv,s ...
  ,Xint ...
  ,Axint,Bxint,Cxint ...
  ,Ayint,Byint,Cyint ...
  ,vlm ...
  ,Istop,Ido,action] ...
  ...
  = prd_2d (NSG,XI,YI,Ich1,thmax,Ich2,spmax,Ich3,spmin);

  if(Istop==1) break; end

%--
end
%--
 
 if(Istop==1) break; end
 
 NSG

 for i=1:NSG+1
  XIT(step,i) = XI(i);
  YIT(step,i) = YI(i);
 end

 NSEG(step)=NSG;

%---
% drop volume
%---

%vlm=0.0;

%for i=1:NSG
% DX = XI(i+1)-XI(i);
% vlm=vlm + YI(i)*YI(i)*DX;
%end

vlm=-vlm/VD;

volume(step)=vlm;

%--
% IRK loop
%---

%-.-.-.
for JRK=1:IRK
%-.-.-.

%---
% geometry
%---

[vnx,vny,crv,s ...
,Xint ...
,Axint,Bxint,Cxint ...
,Ayint,Byint,Cyint ...
,vlm] = splc_geo (NSG,XI,YI);

%---
% mean curvature
%---

crvm(1)     = -crv(1);
crvm(NSG+1) = -crv(NSG+1);

for i=2:NSG
 crvm(i)=-0.5*(crv(i)-vny(i)/YI(i));
end

%---
% jump in normal traction
%---

for i=1:NSG+1
  Dfn(i) = gamma*2.0*crvm(i);
  Dfn(i) = XI(i)*gac*(rhod-rhoa);
  Dfn(i) = gamma*2.0*crvm(i) + XI(i)*gac*(rhod-rhoa);
end

%---
% plot
%---

Iskip=0;
Iskip=1;

if(Iskip==0)
 figure(3)
 plot(s,crv)
 hold on
 plot(s,vnx)
 plot(s,vny)
 plot(s,crvm,'o-')
 plot(s,Dfn,'r')
end

%==============================

%---
% single layer over the interface
% computed over the interface
%---

cf= -1.0/(8.0*pi);

%---
 for j=1:NSG+1  % run over nodes
%---

 X0   =  XI(j);
 Y0   =  YI(j);
 Dfn0 = Dfn(j);

 sum1 = 0.0;  % to test identities
 sum2 = 0.0;

 vx = 0.0;
 vy = 0.0;

%---
 for i=1:NSG
%---

 X1 = XI(i);
 Y1 = YI(i);
 X2 = XI(i+1);
 Y2 = YI(i+1);

 Ising=0;

%---
 if(Ichoose==1)  % lines
%---

   if(j==i | j==i+1)
     Ising=1;
   end

 [Qxx,Qxs ...
 ,Qsx,Qss] = drop_ax_slp_line  ...
 ...
  (X0,Y0 ...
  ,X1,Y1 ...
  ,X2,Y2 ...
  ,NGL ...
  ,Iflow ...
  ,wall ...
  ,sc      ...
  ,RL      ...
  ,Nsum,Np ...
  ,Ising);

  vnx0=0.5*(vnx(i)+vnx(i+1));
  vny0=0.5*(vny(i)+vny(i+1));
  vnx0=-YI(i+1)+YI(i);
  vny0= XI(i+1)-XI(i);
  norm = sqrt(vnx0^2+vny0^2);
  vnx0 = vnx0/norm;
  vny0 = vny0/norm;

  sum1 = sum1 + Qxx*vnx0+Qxs*vny0;
  sum2 = sum2 + Qsx*vnx0+Qss*vny0;

  Dfx0=0.5*(Dfn(i)*vnx(i)+Dfn(i+1)*vnx(i+1));
  Dfy0=0.5*(Dfn(i)*vny(i)+Dfn(i+1)*vny(i+1));

  vx = vx + Qxx*Dfx0 + Qxs*Dfy0;
  vy = vy + Qsx*Dfx0 + Qss*Dfy0;

%---
 else  % spline
%---

  Ising =1;

  Ax=Axint(i); Bx=Bxint(i); Cx=Cxint(i);
  Ay=Ayint(i); By=Byint(i); Cy=Cyint(i);

  Xint1 = Xint(i); Xint2 =Xint(i+1);
  Dfn1  =  Dfn(i); Dfn2  = Dfn(i+1);
  vnx1  =  vnx(i); vnx2   =vnx(i+1);
  vny1  =  vny(i); vny2   =vny(i+1);

  [Qx,Qy ...
  ,Wx,Wy] = drop_ax_slp_spline ...
  ...
  (X0,Y0 ...
  ,X1,Y1 ...
  ,X2,Y2 ...
  ,NGL ...
  ,Xint1,Xint2 ...
  ,Ax,Bx,Cx ...
  ,Ay,By,Cy ...
  ,Dfn1,Dfn2 ...
  ,vnx1,vnx2 ...
  ,vny1,vny2 ...
  ,Dfn0 ...
  ,Iflow ...
  ,wall ...
  ,sc      ...
  ,RL      ...
  ,Nsum,Np ...
  ,Ising ...
  );

  sum1 = sum1 + Qx;
  sum2 = sum2 + Qy;

  vx = vx + Wx;
  vy = vy + Wy;

%---
 end
%---

%---
  end % of interfacial segments
%---

  sum1=sum1*cf;  % should be zero
  sum2=sum2*cf;

%  [sum1 su%2]
%  pause

  vx = vx*cf;
  vy = vy*cf;

  velx(j) = vx;
  vely(j) = vy;

%---
end
%---

vely(1)    =0.0;
vely(NSG+1)=0.0;

%[velx vely]'
%pause

%[vnx vny]'
%pause


if(Iflow==1)
 velx=velx-terminal;
end

%============================

%---
% plot
%---

Ido=0;

if(Ido==1)
 figure(4)
 hold on
 plot(XI(1:NSG+1),velx)
 plot(XI(1:NSG+1),vely,'r')
end

 for i=1:NSG+1
   veln(i) = vnx(i)*velx(i)+vny(i)*vely(i);
 end

 if(JRK==1 & IRK==2)
  XIsave=XI;
  YIsave=YI;
  velxsave=velx;
  velysave=vely;
  velnsave=veln;
  vnxsave=vnx;
  vnysave=vny;
 end

%---
 if(JRK==1)
%---

 for i=1:NSG+1
  if(move==0)
    XI(i) = XI(i)+Dt*velx(i);
    YI(i) = YI(i)+Dt*vely(i);
  else
    XI(i) = XI(i)+Dt*veln(i)*vnx(i);
    YI(i) = YI(i)+Dt*veln(i)*vny(i);
  end
 end

%---
 else
%---

 for i=1:NSG+1
  if(move==0)
    XI(i) = XIsave(i)+0.5*Dt*(velxsave(i)+velx(i));
    YI(i) = YIsave(i)+0.5*Dt*(velysave(i)+vely(i));
  else
    XI(i) = XIsave(i)+0.5*Dt*(velnsave(i)*vnxsave(i)+veln(i)*vnx(i));
    YI(i) = YIsave(i)+0.5*Dt*(velnsave(i)*vnysave(i)+veln(i)*vny(i));
  end
 end

%---
 end
%---

%---
% print
%---

fprintf(file1,'%4d',NSG);
fprintf(file1,'\n');
fprintf(file1,'%5d %12.8f %12.8f %12.8f %12.8f' ...
      ,step,time(step),Dt,volume(step));
fprintf(file1,'\n');
fprintf(file1,'\n');
for i=1:NSG+1
 fprintf(file1,'%13.8f %13.8f %13.8f %13.8f %13.8f %13.8f' ...
 ,XI(i),YI(i),crv(i),velx(i),vely(i));
 fprintf(file1,'\n');
end
fprintf(file1,'\n');

%---
% one more step
%---

time(step+1)=time(step)+Dt;

%-.-.-.
end % of JRK
%-.-.-.

%=========
end
%=========

fprintf(file1,'%5d',999);
fclose(file1);

%---
% save the last profile
%---

 for i=1:NSG+1
  XIT(Nstep+1,i) = XI(i);
  YIT(Nstep+1,i) = YI(i);
 end
 NSEG(Nstep+1)=NSG;

%---
% plot
%---

Ido =0;

if(Ido==1)

figure(99)
hold on
axis equal

many=Nstep+1;

for step=1:many
 if(step==1)
  ntp = NSEG(step)+1;
  plot(XIT(step,1:ntp), YIT(step,1:ntp),'r.-')
  plot(XIT(step,1:ntp),-YIT(step,1:ntp),'r.-')
 elseif(step==many)
  ntp = NSEG(step)+1;
  plot(XIT(step,1:ntp), YIT(step,1:ntp),'c.-')
  plot(XIT(step,1:ntp),-YIT(step,1:ntp),'c.-')
 else
  ntp = NSEG(step)+1;
  plot(XIT(step,1:ntp), YIT(step,1:ntp),'b.-')
  plot(XIT(step,1:ntp),-YIT(step,1:ntp),'b.-')
 end
end

end
