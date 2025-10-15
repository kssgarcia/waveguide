close all
clear all

%=========================
% axisymmetric flow past or due to the
% translation of a particle 
% in free space, in a semi-infinite domain
% or inside a circular tube
%=========================

%=====
Igeom=2;  % ellipsoid
Igeom=3;  % cylinder
Igeom=1;  % sphere
%=====

%=====
Iflow = 2;  % semi-infinite flow
Iflow = 1;  % free space
Iflow = 3;  % inside a cylinder of radius sc
%=====

%---
% initialize
%---

 NGL   = 6;
 mu    = 1.0;
 Vx    = 1.0;
 wall  = 0.0;

 xcnt  = 0.0;
 ycnt  = 0.0;
 sc    = 1.0;
 RL    = 1.4;

 Nsum  = 2;
 Np    = 2;
 rad   = 0.5;
 rad   = 0.4;

 amaj  = 0.0;
 amin  = 0.0;
 tilt  = 0.0;

%=======
% sphere
%=======

if(Igeom==1)

 nelm = 64;
 nelm = 16;
 nelm = 32;

 Itype = 2;   % circular arcs

%---
% element end points
%---

 ratio  = 0.1;
 angle1 = 0.0;
 angle2 = pi;
 sinit  = 0.0;
 Isym   = 1;

 [xe,ye,te,se,xm,ym,tm,sm]...
  ...
  = elm_arc (nelm,ratio,xcnt,ycnt ...
            ,rad,angle1,angle2 ...
            ,sinit,Isym);
%---
% element points
%---

 for i=1:nelm+1
   tel(i) = te(i);
   xel(i) = xe(i);
   yel(i) = ye(i);
   sel(i) = se(i);
 end

%---
% collocation points
%---

 for i=1:nelm
   tcl(i) = tm(i);
   xcl(i) = xm(i);
   ycl(i) = ym(i);
   scl(i) = sm(i);
 end

%---
% element arc length
%---

 for i=1:nelm
  lel(i) = rad*(tel(i+1)-tel(i));
 end

% a = rad;

end

%=======
% spheroid
%=======

if(Igeom==2)

 nelm = 64;
 nelm = 32;
 nelm = 16;

 a = 1.0;
 xcnt = 3.0*a;
 ycnt = 0.0;

 amaj  = a;   % major axis
 amin  = a;   % minor axis
 tilt  = 0.0; % tilting angle

 Itype = 1;   % straight elements
 Itype = 3;   % elliptical segments

%---
% element end points
%---

 Dtheta = pi/nelm; 

 for i=1:nelm+1
   tel(i) = (i-1.0)*Dtheta;
   xel(i) = xcnt+amaj*cos(tel(i));
   yel(i) = ycnt+amin*sin(tel(i));
 end

%---
% collocation points
%---

 for i=1:nelm
  tcl(i) = (i-0.5)*Dtheta;
  xcl(i) = xcnt+amaj*cos(tcl(i));
  ycl(i) = ycnt+amin*sin(tcl(i));
 end

%---
% element arc length
%---

 for i=1:nelm
  lel(i) = sqrt((xel(i+1)-xel(i))^2+(yel(i+1)-yel(i))^2);
 end

%---
end
%---

%=======
% cylinder
%=======

if(Igeom==3)

 Itype = 1;   % straight elements

 a = 0.5;
 b = 0.05;
 b = 0.1;

 nelm1 = 16;
 nelm2 = 64;
 nelm3 = 16;

 Ic = 0;

 x1=a; y1=0; x2=a; y2=b;
 sinit = 0.0;
 ratio1 = 0.1;
 Isym = 0;

 [xe,ye,se,xm,ym,sm]...
 ...
   = elm_line (nelm1,ratio1,x1,y1,x2,y2,sinit,Isym);

 for i=1:nelm1
   Ic = Ic+1;
   tel(Ic) = 0.0;
   xel(Ic) = xe(i);
   yel(Ic) = ye(i);
   sel(Ic) = se(i);
 end

 sinit = se(nelm1+1);
 x1=a; y1=b; x2=-a; y2= b;
 ratio2 = 10;
 Isym = 1;

 [xe,ye,se,xm,ym,sm]...
 ...
   = elm_line (nelm2,ratio2,x1,y1,x2,y2,sinit,Isym);

 for i=1:nelm2
   Ic = Ic+1;
   tel(Ic) = 0.0;
   xel(Ic) = xe(i);
   yel(Ic) = ye(i);
   sel(Ic) = se(i);
 end

 sinit = se(nelm2+1);
 x1 = -a; y1 = b; x2 = -a; y2 = 0;
 ratio3 = 10.0;
 Isym = 0;

 [xe,ye,se,xm,ym,sm]...
 ...
   = elm_line (nelm3,ratio3,x1,y1,x2,y2,sinit,Isym);

 for i=1:nelm3+1
   Ic = Ic+1;
   tel(Ic) = 0.0;
   xel(Ic) = xe(i);
   yel(Ic) = ye(i);
   sel(Ic) = se(i);
  end

  nelm = nelm1+nelm2+nelm3;

 for i=1:nelm
  scl(i) = 0.5*(sel(i)+sel(i+1));
 end

%---
% collocation points
%---

 for i=1:nelm
  tcl(i) = 0.0;
  xcl(i) = 0.5*(xel(i+1)+xel(i));
  ycl(i) = 0.5*(yel(i+1)+yel(i));
 end

%---
% element length
%---

 for i=1:nelm
  lel(i) = sqrt((xel(i+1)-xel(i))^2+(yel(i+1)-yel(i))^2);
 end

%---
end
%---

%---
% plot the elements
%---

 figure(1)
 hold on
 axis equal

 plot(xel, yel,'-k.')
 plot(xel,-yel,'-k.')
 plot(xcl,ycl,'gs')

 if(Iflow==2)
  plot([0 4*a],[0,0],'k--')
  plot([0 0],[-2*a 2*a],'k--')
  axis([-0.5*a 5*a -2*a 2*a])
  plot([wall wall],[-2*a 2*a],'ro')
 end

 if(Iflow==3)
  plot(xel+RL, yel,'-k')
  plot(xel+RL,-yel,'-k')
  plot(xel-RL, yel,'-k')
  plot(xel-RL,-yel,'-k')
  plot([-2 2],[sc, sc],'k-')
  plot([-2 2],[-sc, -sc],'k-')
  plot([-2 2],[0, 0],'k--')
 end

 xlabel('x','fontsize',15)
 ylabel('y','fontsize',15)
 set(gca,'fontsize',15)
 box on

%---
% influence matrix
%---

for i=1:nelm

   X0 = xcl(i);
   Y0 = ycl(i);
   T0 = tcl(i);

   for j=1:nelm

     X1 = xel(j);
     Y1 = yel(j);
     T1 = tel(j);

     X2 = xel(j+1);
     Y2 = yel(j+1);
     T2 = tel(j+1);

     Ising = 0;
     if(i==j) Ising = 1; end

      [Qxx,Qxy ...
      ,Qyx,Qyy] = ...
  ...
      prtcl_ax_slp ...
  ...
        (Iflow ...
        ,X0,Y0,T0 ...
        ,X1,Y1,T1 ...
        ,X2,Y2,T2 ...
        ,NGL ...
        ,Ising ...
        ,Itype ...
        ,xcnt,ycnt ...
        ,rad ...
        ,amaj,amin,tilt ...
        ,wall ...
        ,sc      ...
        ,RL      ...
        ,Nsum,Np ...
        );

    MAT(i,     j)      = Qxx;
    MAT(i,     j+nelm) = Qxy;
    MAT(i+nelm,j)      = Qyx;
    MAT(i+nelm,j+nelm) = Qyy;

    end

    RHS(i) = 8.0*pi*mu*Vx;
    RHS(i+nelm) = 0.0;

   end

  for i=1:nelm
   test = 0.0;
   for j=1:nelm
    test = test+MAT(i,j);
   end
%  test
%  pause
  end

%---
% solve the linear system
%---

  sln = RHS/MAT';

  for i=1:nelm
    fclx(i) = sln(i); 
    fcly(i) = sln(i+nelm);
  end

%---
% force on the particle
%---

  forcex = 0.0;

  for i=1:nelm
   forcex = forcex + fclx(i)*2*pi*ycl(i)*lel(i);
  end

  if(Igeom==1)
   force_coefficient = forcex/(6.0*pi*mu*rad)
  elseif(Igeom==2)
   force_coefficient = forcex/(6.0*pi*mu*a)
  elseif(Igeom==3)
   force_coefficient = forcex/(6.0*pi*mu*a)
  end

%---
% pressure drop in a tube
%---

 if(Iflow==3)

 pd = 0.0;

  for i=1:nelm
   pd = pd + fclx(i)*2*pi*(sc^2-ycl(i)^2)*ycl(i)*lel(i);
  end
 
  pd = 2.0*pd/pi

 end


%---
% plot the solution
%---

if(Igeom==1)
   figure(2)
   hold on
   plot(xcl,fclx,'ko-')
   plot(xcl,fcly,'ko--')
   xlabel('x','fontsize',15)
   ylabel('f_x, f_y','fontsize',15)
   set(gca,'fontsize',15)
end

if(Igeom==2)
   figure(2)
   hold on
   plot(xcl,fclx,'ko-')
   plot(xcl,fcly,'ko--')
   xlabel('x','fontsize',15)
   ylabel('f_x, f_y','fontsize',15)
   set(gca,'fontsize',15)
end

if(Igeom==3)

   alpha = 2*pi/(log(2*a/b)-0.5);
   slender = alpha/(2*pi*b);

   figure(2)
   hold on
   set(gca,'fontsize',15)
   plot(scl,fclx,'k-')
   plot([b,b+2*a],[slender,slender],'k--')
   xlabel('s/a','fontsize',15)
   ylabel('f_x','fontsize',15)
   axis([0,2*(a+b),0, 20])
   box on

   figure(3)
   hold on
   set(gca,'fontsize',15)
   plot(scl,fcly,'k-')
   xlabel('s/a','fontsize',15)
   ylabel('f_\sigma','fontsize',15)
   axis([0,2*(a+b),-20, 20])
   box on

end
