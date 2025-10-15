close all
clear all

%++++++++++++++++++++++++++
% Numerical simulation of the
% two-dimensional capsule deformation
% in shear flow
%++++++++++++++++++++++++++

%---
% globals
%---

global NGL
global Iflow
global wall 
global period
global shrt
global mu1
global mu2

%---
Igraph1 = 0;
%---

%======
% parameters
%======

a = 1.0; % equivalent radius

%---
% viscosities
%---

mu1 = 1.0;

mu2 = 100.0;
mu2 = 20.0;
mu2 = 1.0;
mu2 = 5.0;

%---
% choose the flow
%---

Iflow = 3;  % periodic above a wall
Iflow = 0;  % no flow
Iflow = 1;  % infinite shear flow
Iflow = 2;  % solitary above a wall

wall = 0.0;   % wall at y=wall
period = 3.5; % for Iflow=3

shrt = 0.0;  % shear rate
shrt = 1.0;  % shear rate

if(Iflow==0)
 shrt = 0.0;
end

gamma = 0.0;  % surface tension

%---
% elastic modulus
%---

elst = 10.0;
elst = 5.0;
elst = 1.0;
elst = 50.0;
elst = 20.0;

%---
% Hat E_B
%---

bmodhat = 0.010;
bmodhat = 0.000;
bmodhat = 0.001;
bmodhat = 0.050;
bmodhat = 0.005;

bmod = bmodhat*elst*a*a;

%bmod = 0.00;
%bmod = 0.05;  % bending modulus
%bmod = 0.02;  % bending modulus

NSG = 128;
NSG = 96;
NSG = 32;
NSG = 64;

XDC=0.0;  % initial cell center position
YDC=0.0;

if(Iflow==2)
 YDC=3.00;
 YDC=4.00;
 YDC=1.25;
 YDC=3.00;
 YDC=1.50;
 YDC=2.00;
end

if(Iflow==3)
 YDC=1.5;
end

Ishape = 1; % ellipse
Ishape = 2; % rbc

rotangle = 0.0;
rotangle = 3.0*pi/4.0;
rotangle = pi/4.0;
rotangle = pi/8.0;

NGL = 2;  % Gauss-Legendre
NGL = 6;  % Gauss-Legendre

IRK = 1;
IRK = 2;

Nsteps = 400;
Nsteps = 500;
Nsteps = 300;
Nsteps = 1;
Nsteps = 2000;
Nsteps = 600;
Nsteps = 1000;

%---
tfinal = 8.0;
tfinal = 12.0;
tfinal = 25.0;
tfinal = 200.0;
%---

Dt = 0.015;
Dt = 0.01;
Dt = 0.02;

%---
% smoothing the tensions
%---

Nsmoothdf = 10;

%---
% point redistribution
%---

Ich1=1;
Ich2=1;
Ich3=1;

%---
thmax=0.25*pi;

%---
spmax=0.10;
spmin=0.05;
%---

%---
spmax=0.15;
spmin=0.075;

%---
spmax=0.20;
spmin=0.10;
%---


%---
spmax=0.15;
spmin=0.05;
%---

%~~~~~~~~~~~~
% end of input
%~~~~~~~~~~~~

%---
% prepare
%---

Dth = 2*pi/NSG;
csrot = cos(rotangle);
snrot = sin(rotangle);

%---
% shape: ellipse
%---

if(Ishape==1)

 axis1 = 1.0;
 axis2 = 2.0;

 for i=1:NSG+1
  th = Dth*(i-1);
  X(i) = axis1*a*cos(th);
  Y(i) = axis2*a*sin(th);
 end

 scale = 1.0;

end

%---
% shape: RBC
%---

if(Ishape==2)

  for i=1:NSG+1
%    tt = tt+3.00D0*Dth*sin(2.0D0*tt);
    tt = (i-1.0D0)*Dth;
    cs = cos(tt);
    sn = sin(tt);
    X(i) = cs;
    Y(i) = 0.5*sn*(0.207+2.003*cs^2-1.123*cs^4);
  end

  scale = 1.3858189*a;

end

%---
% rotate, expand, and translate
%---

for i=1:NSG+1
 Xsave = X(i);
 Ysave = Y(i);
 X(i) = Xsave*csrot - Ysave*snrot;
 Y(i) = Xsave*snrot + Ysave*csrot;
 X(i) = XDC + X(i)*scale;
 Y(i) = YDC + Y(i)*scale;
end

%hold on
%plot(X,Y,'-')
%plot(X,Y,'ro')
%axis equal

%---
% unstressed geometry
%---

 [vnx,vny,crvuns,suns ...
 ,Xint ...
 ,Axint,Bxint,Cxint ...
 ,Ayint,Byint,Cyint ...
 ,areauns,centerx,centery ...
 ,aspectuns,angle1,angle2] = splc_geo (NSG,X,Y);

%figure
%plot(s,vnx) 
%plot(s,crv,'-o') 

%---
% prepare to run
%---

time = 0.0;
file1 = fopen('rbc_2d.out','w');

%---
% continue a calculation?
%---

Icontinue = 1;  % yes
Icontinue = 0;

if(Icontinue==1)

 clear X Y crvuns suns

 [NSG,X,Y,suns,crvuns,time] = rbc_input();

 X = smooth(NSG,X,Nsmooth);
 Y = smooth(NSG,Y,Nsmooth);

 figure
 hold on
 plot(X(1:NSG+1),Y(1:NSG+1),'ro')
 pause

end

%---
% a circle (for relaxation)
%---

Icircle = 1;
Icircle = 0;

if(Icircle==1)

 radius = sqrt(areauns/pi);

 for i=1:NSG+1
  th = Dth*(i-1);
  X(i) = radius*cos(th);
  Y(i) = radius*sin(th);
 end

 shrt = 0.0;
 scale = 1.0;
 elst = 1.0;

 bmodhat = 0.005;
 bmodhat = 0.000;

 bmod = bmodhat*elst*a*a;

 Dt = 0.5;

 Nsteps = 200;
 tfinal = 99999.0;

 if(bmodhat<0.0000001)
   Nsteps = 50;
%  Ich1=0;Ich2=0;Ich3=0;
 end

end

%==================
% time stepping
%==================

for step=1:Nsteps+1

%---
% smoothing
%---

%Nsmooth = 1;
%Nsmooth = 0;
%if(step==100 | step==200 | step==300 | step==400 ...
%  |step==500 | step==600 | step==700 | step==800 ...
%  |step==900 | step==1000)
%  Nsmooth = 10;
%  Nsmooth = 0;
%end

Nsmooth = 2;
Nsmooth = 4;
Nsmooth = 0;

if(Nsmooth>0)
 X = smooth(NSG,X,Nsmooth);
 Y = smooth(NSG,Y,Nsmooth);
end

%---
% point redistribution
%---

Isym = 1;
Isym = 0;

Irepeat = 1;

while(Irepeat==1)

 [NSG,X,Y ...
 ,vnx,vny ...
 ,crv,s,Xint ...
 ,Axint,Bxint,Cxint ...
 ,Ayint,Byint,Cyint ...
 ,area,centerx,centery ...
 ,suns,crvuns ...
 ,Istop,Irepeat,action] ...
 ...
 = prd_2d (NSG,X,Y ...
          ,suns,crvuns ...
          ,Ich1,thmax ...
          ,Ich2,spmax ...
          ,Ich3,spmin,Isym);

if(Istop==1) break; end

end

if(Istop==1) break; end

%---
% normalize the area
%---

 Inormalize=1;
 Inormalize=0;

 if(Inormalize==1)
  factor = sqrt(areauns/area);
  for i=1:NSG+1
    X(i) = (X(i)-centerx)*factor+centerx;
    Y(i) = (Y(i)-centery)*factor+centery;
  end
 end

%---
% geometry
%---

 [vnx,vny,crv,s ...
 ,Xint ...
 ,Axint,Bxint,Cxint ...
 ,Ayint,Byint,Cyint ...
 ,area,centerx,centery ...
 ,aspect,angle1,angle2] = splc_geo (NSG,X,Y);


Ido=1;
Ido=0;

if(Ido==1)
 figure
 plot(X(1:NSG+1),Y(1:NSG+1),'o')
 hold on
 plot(X(1),Y(1),'+')
 axis equal
 figure
 hold on
 plot(s(1:NSG+1),crv(1:NSG+1),'+-')
 plot(s(1:NSG+1),crvuns(1:NSG+1),'ro-')
 figure
 hold on
 plot(s(1:NSG+1),suns(1:NSG+1),'ro-')
 pause
end

%---
% animation
%---

if(step==1)
  Handle100 = plot(X(1:NSG+1),Y(1:NSG+1),'k.','LineWidth',2);
  hold on
  Handle2 = plot(X(1:NSG+1),Y(1:NSG+1),'ro');
% Handle9 = patch(X(1:NSG+1),Y(1:NSG+1),'r');
  if(Iflow==3)
   Handle3 = plot(X(1:NSG+1),Y(1:NSG+1),'ro');
   Handle4 = plot(X,Y,'b-.');
   set(Handle3,'EraseMode','xor')
   set(Handle4,'EraseMode','xor')
  end
  set(Handle100,'EraseMode','xor')
   set(Handle2,'EraseMode','xor')
%  set(Handle9,'EraseMode','xor')
  xlabel('x','fontsize',15)
  ylabel('y','fontsize',15)
  set(gca,'fontsize',15)
  if(Iflow==1)
    axis([-2.0 2.0 -2.0 2.0])
  end
  axis equal
  axis off
end

 set(Handle100,'XData',X(1:NSG+1),'YData',Y(1:NSG+1))
%set(Handle2,'XData',X(1:NSG+1),'YData',Y(1:NSG+1));
%set(Handle9,'XData',X(1:NSG+1),'YData',Y(1:NSG+1));
  if(Iflow==3)
   set(Handle3,'XData',X+period,'YData',Y)
   set(Handle3,'EraseMode','xor')
   set(Handle4,'XData',X+period,'YData',Y)
   set(Handle4,'EraseMode','xor')
  end
 drawnow
 pause(0.1)

 if(Iflow==2)
  displ = 0.8*step*Dt*shrt*YDC;
%  axis([-2+displ 2+displ -0.5 3.5])
%  plot([-2+displ 2+displ],[0 0],'LineWidth',2);
 end

 if(Iflow==3)
  displ = 0.8*step*Dt*shrt*YDC;
  axis([-2+displ 2+displ -0.5 3.5])
  plot([-2+displ 2+displ],[0 0],'LineWidth',2);
 end

%---
% tension
%---

 for i=1:NSG+1
  srtn(i) = gamma;
 end

%---
% Df
%---

 [Dfx,Dfy,elten,tsten] = rbc_df ...
...
    (NSG ...
    ,Xint ...
    ,srtn ...
    ,suns,crvuns ...
    ,s,crv,vnx,vny ...
    ,elst ...
    ,bmod ...
    ,Nsmoothdf);

if(Igraph1==1)
 figure(5)
 hold on
 plot(s(1:NSG+1),Dfx,'-') 
 figure(6)
 hold on
 plot(s(1:NSG+1),Dfy,'r--') 
end

%---
% velocity
%---

 [U,V] = rbc_velocity ...
 ...
  (NSG,X,Y,Xint ...
  ,Axint,Bxint,Cxint ...
  ,Ayint,Byint,Cyint ...
  ,Dfx,Dfy ...
  ,vnx,vny);

%---
% tangential velocity
%---

 for i=1:NSG+1
  tanvel(i) = U(i)*vny(i)-V(i)*vnx(i);
 end

%---
% record
%---

 angle1=angle1/pi
 angle2=angle2/pi
 aspect=aspect/aspectuns

 fprintf(file1,'%4d',NSG);
 fprintf(file1,'\n');
 fprintf(file1,'%5d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f',step,time,Dt,area,centerx,centery,aspect,angle1,angle2);
 fprintf(file1,'\n');
 for i=1:NSG+1
  fprintf(file1,'%3d %13.8f %13.8f %13.8f %13.8f %13.8f %13.8f %13.8f %13.8f %13.8f %13.8f %13.8f %13.8f %13.8f %13.8.f' ...
 ,i,X(i),Y(i),s(i),crv(i),U(i),V(i) ...
   ,tanvel(i),elten(i),tsten(i),suns(i),crvuns(i));
  fprintf(file1,'\n');
 end
 fprintf(file1,'\n');

  if(time>tfinal)
   break
  end

  if(step>Nsteps+1)
   break
  end

%---
% advance
%---

 if(IRK==2)
  for i=1:NSG+1
   Xsave(i) = X(i);
   Ysave(i) = Y(i);
   Usave(i) = U(i);
   Vsave(i) = V(i);
  end
 end

 for i=1:NSG+1
  X(i) = X(i)+Dt*U(i);
  Y(i) = Y(i)+Dt*V(i);
 end

 time = time+Dt;
 time

%---
if(IRK==2)
%---

 [vnx,vny,crv,s ...
 ,Xint ...
 ,Axint,Bxint,Cxint ...
 ,Ayint,Byint,Cyint ...
 ,area,centerx,centery ...
 ,aspect,angle1,angle2] = splc_geo (NSG,X,Y);

 [Dfx,Dfy,elten,tsten] = rbc_df ...
...
    (NSG ...
    ,Xint ...
    ,srtn ...
    ,suns,crvuns ...
    ,s,crv,vnx,vny ...
    ,elst ...
    ,bmod ...
    ,Nsmoothdf);

 [U,V] = rbc_velocity ...
 ...
  (NSG,X,Y,Xint ...
  ,Axint,Bxint,Cxint ...
  ,Ayint,Byint,Cyint ...
  ,Dfx,Dfy ...
  ,vnx,vny);

 for i=1:NSG+1
  X(i) = Xsave(i)+0.5*Dt*(U(i)+Usave(i));
  Y(i) = Ysave(i)+0.5*Dt*(V(i)+Vsave(i));
 end

%---
end  % of IRK2
%---

%=================
end % of time steps
%=================

fprintf(file1,'999');
fprintf(file1,'\n');
fprintf(file1,'\n');
fprintf(file1,'mu1=%12.8f;',mu1);
fprintf(file1,'\n');
fprintf(file1,'mu2=%12.8f;',mu2);
fprintf(file1,'\n');
fprintf(file1,'elst=%12.8f;',elst);
fprintf(file1,'\n');
fprintf(file1,'gamma=%12.8f;',gamma);
fprintf(file1,'\n');
fprintf(file1,'wall=%12.8f;',wall);
fprintf(file1,'\n');
fprintf(file1,'a=%12.8f;',a);
fprintf(file1,'\n');
fprintf(file1,'bmodhat=%12.8f;',bmodhat);
fprintf(file1,'\n');
fprintf(file1,'bmod=%12.8f;',bmod);
fprintf(file1,'\n');
fprintf(file1,'XDC=%12.8f;',XDC);
fprintf(file1,'\n');
fprintf(file1,'YDC=%12.8f;',YDC);
fprintf(file1,'\n');
fprintf(file1,'shrt=%5.0f;',shrt);
fprintf(file1,'\n');
fprintf(file1,'NSG=%4d;',NSG);
fprintf(file1,'\n');
fprintf(file1,'IRK=%4d;',IRK);
fprintf(file1,'\n');
fprintf(file1,'Dt=%12.8f;',Dt);
fprintf(file1,'\n');
fprintf(file1,'Iflow=%4d;',Iflow);
fprintf(file1,'\n');
fprintf(file1,'Ishape=%4d;',Ishape);
fprintf(file1,'\n');
fprintf(file1,'NGL=%4d;',NGL);
fprintf(file1,'\n');
fprintf(file1,'period=%5.0f;',period);
fprintf(file1,'\n');
fprintf(file1,'thmax=%12.8f;',thmax);
fprintf(file1,'\n');
fprintf(file1,'spmax=%12.8f;',spmax);
fprintf(file1,'\n');
fprintf(file1,'spmin=%12.8f;',spmin);
fprintf(file1,'\n');
fprintf(file1,'Inormalize=%4d;',Inormalize);
fprintf(file1,'\n');
fprintf(file1,'Nsmooth=%4d;',Nsmooth);
fprintf(file1,'\n');
fprintf(file1,'Nsmoothdf=%4d;',Nsmoothdf);
%---

fclose(file1);
