close all
clear all

%============================================
% FDLIB, BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved
%
% This program is to be used only under the
% stipulations of the licensing agreement
%============================================

%--------------------------------------------------
% Dynamical simulation of the sloshing of an inviscid
% liquid inside a (two-dimensional) rectangular tank
%
%  LEGEND:
%  ------
%
%  N:   Number of segments along the free surface
%
%  x ,y :  coordinates of points along the free surface
%  xs,ys:  saved coordinates
%  phi:    potential at nodes
%  ps:     saved potential at nodes
%  crv:    curvature at nodes
%
%  wall1: One vertical wall located at x = wall1
%  wall2: One vertical wall located at x = wall2
%  wall3:   Horizontal wall located at y = wall3
%
%  NGL:   Number of Gauss--Legendre points for numerical integration
%         over the free-surface elements
%
%  xm,ym:    Coordinates of mid-points
%  vnx,vny:  Coordinates of normal vector at mid-points
%
%  rho:   fluid density
%  gamma: surface tension
%  rmu:   viscous damping coefficient
%
%  h:     unperturbed free surface located at y = h
%  a0:    initial amplitude of the free surface
%  a0phi: initial amplitude of the potential along the free surface
%
%  accx: x-component of the container acceleration
%  accy: y-component of the container acceleration
%
%  accxtime: duration of the x-component of the container acceleration
%  accytime: duration of the y-component of the container acceleration

%  Unm:  normal velocity at mid-points
%  phim: potential at mid-points
%
%  slp(i,j):  integral of the slp at ith point over jth segment
%  dlp(i,j):  integral of the dlp at ith point over jth segment
%
%  al:  arc-length of segments
%
%  AB * SLN = BM: Linear system for the normal velocity
%
%  Nprint:   Will print a profile after Nprint time steps
%  Nsmooth:  Will smooth profile and potential after Nsmooth steps
%
%  gac:   acceleration of gravity
%
%  Move = 0    marker points move with total velocity
%         1    marker points move with normal velocity
%
%  Iread = 0   Will generate the initial condition
%          1   Will read from file: tank_2d.inp
%
%  Iflp: index for the extrapolation of the contact line
%        1 for linear extrapolation
%        2 for quadratic extrapolation
%-------------------------------------------

 global vnxm vnym elml
 global wall1 wall2 wall3
 global NGL ZZ WW

%----------
% constants
%----------

 pi2 = 2.0*pi;
 tot  = 3.0/2.0;

%----------
% input data
%----------

wall1 = -0.5;
wall2 = 0.5;
wall3 = 0.0;
h0 = 1.0;      % initial liquid height
a0 = 0.01;     % disturbance amplitude for free surface
a0phi = 0.00;  % disturbance amplitude for free surface potentian
gac = 1.00;    % gravity acceleration
rho = 1.00;    % density
gamma = 0.00;  % surface tension
rmu = 0.00;    % viscous dissipation coefficient (phenomenological)
N = 16;          % Number of elements
NGL = 6;         % Gauss--Legendre
Move = 1;        % 0 for total velocity, 1 for normal velocity
Iflp = 1;       % index for the positioning of the first and last points
Dt = 0.002;     % time step
accx = 1.0; accy = 0.0 ;   % components of external acceleration
accxtime = 0.4; accytime = 1.0;  % duration of external acceleration
Nstep = 5000;
Nsmooth = 1;    %   After how many steps smoothing
Ired = 1;  % point redistribution flag
Ich1 = 1; thmax = pi/4; %   point redistribution
Ich2 = 1; spmax = 0.10; %   point redistribution
Ich3 = 1; spmin = 0.02; %   point redistribution

%---
% prepare
%---

 width = wall2-wall1;

 N1 = N+1;

 [ZZ,WW] = gauss_leg(NGL); % Gauss--Legendre base points and weights

%-----------------------------
% generate the initial position of the free surface
% distribution of potential (phi)
% over the free surface
%------------------------------

 wn = pi2/width;    %  wave number

 Dx = width/N;

 for i=1:N1
    tmp = (i-1.0D0)*Dx;
    x(i) = wall1+tmp;
    y(i) = wall3 + h0 + a0*cos(wn*tmp);
    phi(i) =         a0phi*cos(wn*tmp);
  end

%-----------
% initialize
%-----------

  time = 0.0;

  Istep   = 1;      % step counter
  Ismooth = 1;      % steps for smoothing

%-----------------
% prepare graphics
%-----------------

  figure(1)
  figure(2)

%-.-.-.-.-.-.-.-.-.-.-.-
% SIMULATION BEGINS HERE
%-.-.-.-.-.-.-.-.-.-.-.-

%-------------
 for istep=1:Nstep
%-------------

 %=========
 if(Ired==1) % point redistribution module
 %=========

 Isym = 0;
 Italk = 1;
 Irepeat = 1;

  while(Irepeat==1)

         [N,x,y ...
         ,vnx,vny ...
         ,crv,s,Xint ...
         ,Axint,Bxint,Cxint ...
         ,Ayint,Byint,Cyint ...
         ,area,centerx,centery ...
         ,phi,phi ...
         ,Istop ...
         ,Irepeat ...
         ,Iaction] ...
...
   = prd_2d_open (N ...
            ,x,y ...
            ,phi,phi ...
            ,Ich1,thmax ...
            ,Ich2,spmax ...
            ,Ich3,spmin ...
            ,Isym ...
            ,Italk);

     if(Istop==1)
      disp(' prd_2d_dr: something went wrong');
      disp('            code:');
      disp(Iaction)
      break
    end

  end

  -area;

  N1 = N+1;

  for i=1:N1
   y(i) = -y(i)*width/area;
  end

 %=========
 end
 %=========

   [istep time]

   if(Ismooth==Nsmooth)
    for smoothloop=1:1

     for i=1:N1        % save old values
       xs(i) =   x(i);
       ys(i) =   y(i);
       ps(i) = phi(i);
     end

     for i=3:N-1
       x(i) = (-xs(i-2)+4.0*xs(i-1)+10.0*xs(i)+4.0*xs(i+1) ...
               -xs(i+2))/16.0;
       y(i) = (-ys(i-2)+4.0*ys(i-1)+10.0*ys(i)+4.0*ys(i+1) ...
               -ys(i+2))/16.0;
       phi(i) = (-ps(i-2)+4.0*ps(i-1)+10.0*ps(i)+4.0*ps(i+1) ...
               -ps(i+2))/16.0;
      end
    end

    Ismooth = 0;

   end

%------
% plot
%-----

  figure(1)
  clf
  hold on
  axis([wall1-0.1,wall2+0.1,wall3-0.1,wall3+2*width])
  axis square
  box on
  plot([wall1,wall1],[0,1.5*h0],'k')
  plot([wall2,wall2],[0,1.5*h0],'k')
  plot([wall1,wall1+width],[wall3,wall3],'k-')
  plot(x(1:N+1),y(1:N+1),'k.-')
  patch([x wall2 wall1 x(1)],[y wall3 wall3 y(1)],'y')
  pause(0.001)

  figure(2)
  clf
  plot(x(1:N+1),phi(1:N+1),'k.-')
  xlabel('x')
  ylabel('\phi')
  pause(0.001)

%------------------------
% Define segment mid-points
%
% Compute: phi at mid-points
%          element length (elml)
%          normal vector pointing into the fluid
%------------------------

  for i=1:N

    xm(i) = 0.5*(x(i+1)+x(i));
    ym(i) = 0.5*(y(i+1)+y(i));
    elml(i) = sqrt((y(i+1)-y(i))^2+(x(i+1)-x(i))^2);

    phim(i) = 0.5*(phi(i+1)+phi(i));

    vnxm(i) =   (y(i+1)-y(i))/elml(i);   % normal vector
    vnym(i) = - (x(i+1)-x(i))/elml(i);   % normal vector

  end

%-----------------------------------------------
% compute the single-layer and double-layer potentials
% at the elements mid-points integrated over the segments
%
% Set up a linear system for the normal velocity
%-----------------------------------------------

      for i=1:N

       x0 = xm(i);
       y0 = ym(i);

       for j=1:N

        j1 = j+1;
        x1 = x(j);
        y1 = y(j);
        x2 = x(j1);
        y2 = y(j1);

        [slp(i,j), dlp(i,j)] = tank_2d_sdlp  ...
       ...
        (x0,y0 ...
        ,i ...
        ,j ...
        ,x1,y1 ...
        ,x2,y2 ...
        );

        AB(i,j) = slp(i,j);      % influence matrix

       end

%---
% compute the right-hand side (rhs) in terms of the dlp
%---

       accum = 0.0;
       for j=1:N
        accum = accum + phim(j)*dlp(i,j);
       end

       BM(i) = accum - 0.5D0*phim(i);     %  right-hand side

      end

%---
% solve the linear system
%---

     sln = BM/AB';
%     sln

%------------------------
% distribute the solution (normal velocity at mid-nodes)
%------------------------

      for i=1:N
        Unm(i) = sln(i);
      end 

%--------------------------------------------
% Variables at element end-nodes
%
% Compute: normal velocity by interpolation (Un)
%          tangential velocity by numerical differentiation
%          normal vector
%          tangential vector 
%          curvature
%--------------------------------------------

  for i=2:N

      dn = elml(i-1)+elml(i);
      Un(i) = (Unm(i-1)*elml(i)+Unm(i)*elml(i-1))/dn;

      x1  = -elml(i-1);
      x2  =  elml(i);
      x21 =  x2-x1;

%---
% quadratic differentiation for phi
% with respect to arc length
%---

      y1  = phi(i-1);
      y0  = phi(i);
      y2  = phi(i+1);
      aa  = ((y2-y0)/x2-(y1-y0)/x1)/x21;
      bb  =  (y2-y0)/x2 - aa*x2;
      Ut(i) = bb;

%---
% quadratic differentiation for x
% with respect to arc length
%---

      y1 = x(i-1);                 % quadratic differentiation
      y0 = x(i);
      y2 = x(i+1);
      aa = ((y2-y0)/x2-(y1-y0)/x1)/x21;
      bb =  (y2-y0)/x2 - aa*x2;
      xp = bb;
      xpp= 2.0D0*aa;

%---
% quadratic differentiation for y
% with respect to arc length
%---

      y1 = y(i-1);           % quadratic differentiation
      y0 = y(i);
      y2 = y(i+1);
      aa = ((y2-y0)/x2-(y1-y0)/x1)/x21;
      bb =  (y2-y0)/x2 - aa*x2;
      yp = bb;
      ypp= 2.0D0*aa;

      tnm    = sqrt(xp*xp+yp*yp);
      tnx(i) = xp/tnm;
      tny(i) = yp/tnm;

      vnx(i) =  tny(i);    % normal vector points into the fluid
      vny(i) = -tnx(i);

%---
% curvature
%---

      crv(i) = (xpp*yp-ypp*xp)/(xp*xp+yp*yp)^tot;

   end

%-------------------------------------------
% time step
%
% update phi using Bernoulli's equation
%-------------------------------------------

  %---
  for i=2:N
  %---

       press  = gamma*crv(i);                % capillary pressure
       Ums    = Un(i)^2 + Ut(i)^2;
       phi(i) = phi(i) - Dt*rmu*phi(i);
       yref   = y(i)-(wall3+h0);

%-----------------------
       if(Move==0)     % total velocity
%-----------------------

        x(i)   = x(i) + Dt*( Un(i)*vnx(i)+Ut(i)*tnx(i) );
        y(i)   = y(i) + Dt*( Un(i)*vny(i)+Ut(i)*tny(i) );
        phi(i) = phi(i) + Dt*(0.5D0*Ums-press/rho-gac*yref);

%----------
       elseif(Move==1)           % normal velocity
%----------

        x(i)   = x(i) + Dt*Un(i)*vnx(i);
        y(i)   = y(i) + Dt*Un(i)*vny(i);
        phi(i) = phi(i) + Dt* (Un(i)*Un(i) ...
                              -0.5D0*Ums ...
                              -press/rho ...
                              -gac*yref);
%------------
       end
%------------

%---
% account for the tank acceleration:
%---

    if(time <= accxtime & accx > 0) 
      phi(i) = phi(i) - Dt*accx*x(i);
    end

    if(time <= accytime & accy > 0)
      phi(i) = phi(i) - Dt*accy*y(i);
    end

  %---
  end   % over points
  %---

  if(time <= accxtime) 
      disp('accelerating in x')
  end
  if(time <= accytime)
      disp('accelerating in y')
  end

  time = time + Dt;

  Ismooth = Ismooth +1;
  Istep   = Istep   +1;

%------------------------------------------
% Compute the y position and phi of first point
%------------------------------------------

   if(Iflp==1)    % linear extrapolation

        x1 = elml(1);
        x2 = x1+elml(2);
        y1 = y(2);
        y2 = y(3);
        y(1) = y1*(0.0-x2)/(x1-x2)+y2*(0.0-x1)/(x2-x1);
        y1 = phi(2);
        y2 = phi(3);
        phi(1) = y1*(0.0-x2)/(x1-x2) ...
                +y2*(0.0-x1)/(x2-x1);

   elseif (Iflp==2)  % quadratic extrapolation

        x1 = elml(1);
        x2 = x1+elml(2);
        x3 = x2+elml(3);
        y1 = y(2);
        y2 = y(3);
        y3 = y(4);
        y(1) = y1*(0.0-x2)*(0.0-x3)/((x1-x2)*(x1-x3)) ...
             + y2*(0.0-x1)*(0.0-x3)/((x2-x1)*(x2-x3)) ...
             + y3*(0.0-x1)*(0.0-x2)/((x3-x1)*(x3-x2));
        y1 = phi(2);
        y2 = phi(3);
        y3 = phi(4);
        phi(1) = y1*(0.0-x2)*(0.0-x3)/((x1-x2)*(x1-x3)) ...
               + y2*(0.0-x1)*(0.0-x3)/((x2-x1)*(x2-x3)) ...
               + y3*(0.0-x1)*(0.0-x2)/((x3-x1)*(x3-x2));

    end

%-----------------------------------------
% Compute the y position and phi of last point
% numbered N+1
%-----------------------------------------

   if(Iflp==1)  % linear extrapolation

        x1 = elml(N);
        x2 = x1+elml(N-1);
        y1 = y(N);
        y2 = y(N-1);
        y(N1) = y1*(0.0-x2)/(x1-x2)+y2*(0.0-x1)/(x2-x1);
        y1 = phi(N);
        y2 = phi(N-1);
        phi(N1) = y1*(0.0-x2)/(x1-x2) ...
                + y2*(0.0-x1)/(x2-x1);

   elseif(Iflp==2)  % quadratic extrapolation

        x1 = elml(N);
        x2 = x1+elml(N-1);
        x3 = x2+elml(N-2);
        y1 = y(N);
        y2 = y(N-1);
        y3 = y(N-2);
        y(N1) = y1*(0.0-x2)*(0.0-x3)/((x1-x2)*(x1-x3)) ...
              + y2*(0.0-x1)*(0.0-x3)/((x2-x1)*(x2-x3)) ...
              + y3*(0.0-x1)*(0.0-x2)/((x3-x1)*(x3-x2));
        y1 = phi(N);
        y2 = phi(N-1);
        y3 = phi(N-2);
        phi(N1) = y1*(0.0-x2)*(0.0-x3)/((x1-x2)*(x1-x3)) ...
                + y2*(0.0-x1)*(0.0-x3)/((x2-x1)*(x2-x3)) ...
                + y3*(0.0-x1)*(0.0-x2)/((x3-x1)*(x3-x2));

    end

  % stop when a point moves outside the container

  Istop = 0;

  for i=2:N
   if(x(i)<wall1) Istop=1; break; end
   if(x(i)>wall2) Istop=1; break; end
   if(y(i)<wall3) Istop=1; break; end
  end

 if(Istop==1) disp('the simulation failed'); break; end

%-------------
 end  % over time steps
%-------------

%---------------------
% simulation has ended
%---------------------

