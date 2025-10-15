function [a,x,s] = drop_ax (Jsp ...
  ,gac,gamma,rhod,rhoa,volume ...
  ,alpha,npts,epsilon,maxiter,tol ...
  )
 
%====================================================
% Hydrostatic shape of an axisymmetric
% sessile drop resting on a horizontal plate,
% or pendant drop hanging underneath a horizontal plate,
% for a specified volume and contact angle.
%====================================================

%--------
% prepare
%--------

 Iflag=0;

 npts1 = npts+1;
 drho  = rhod-rhoa ;            % density difference
 capls = gamma/(gac*abs(drho));  % square of the capillary number

 Isp = 1;                        % Isp is an orientation index
 if(drho<0) 
   Isp = -Isp;
 end
 if(Jsp==-1)
   Isp = -Isp;
 end

%----------------------------------
% To start, assume drop shape
% is a truncated sphere
% and compute the sphere radius "a"
% as a function of 
% volume and contact angle.
%----------------------------------

  cosa = cos(alpha);

  a = (3.0D0*volume/pi)/(2.0D0+cosa^3-3.0D0*cosa);
  a = a^(1.0D0/3.0D0);

  shp(1) = 2.0D0/a;     % twice the mean curvature

%---
% Compute initial solution of the odes
% to start-up the secant method
%---

  dpsi = alpha/npts;

  Ic=1;   % counter

   [x,s,volume_sh] = drop_ax_ode ...
   ...
   (npts ...
   ,capls ...
   ,Isp ...
   ,dpsi ...
   ,shp(Ic) ...
   );

   error(Ic) = volume_sh - volume;
   err = abs(error(Ic));

%-------------------------
% second start-up solution 
%-------------------------

   Ic=2;
   shp(2) = shp(1)+epsilon;

   [x,s,volume_sh] = drop_ax_ode ...
   ...
   (npts ...
   ,capls ...
   ,Isp ...
   ,dpsi ...
   ,shp(Ic) ...
   );

    error(Ic) = volume_sh - volume;
    err = abs(error(Ic));

%---------------------------------------
% iterate on shp using the secant method
% until convergence
%---------------------------------------

   for iter=1:maxiter

    Ic = Ic+1;

%---
% secant updating
%---
     
   Icb = Ic-2;
   Ica = Ic-1;

   dedc = (error(Ica)-error(Icb))/(shp(Ica)-shp(Icb));
   shp(Ic) = shp(Ica)-error(Ica)/dedc;

   [x,s,volume_sh] = drop_ax_ode ...
   ...
   (npts ...
   ,capls ...
   ,Isp ...
   ,dpsi ...
   ,shp(Ic) ...
   );

   error(Ic) = volume_sh - volume;
   err       = abs(error(Ic));

   if(err<tol)
    break
   end

%---
   end
%---

   if(iter==maxiter)
    disp('drop_ax: ODE solver failed')
    Iflag=1;
    return
   end

%---
%  shift to reset the origin at the wall
%---

 shift = x(npts1);
 x=x-shift;

%---
% done
%---
 
 return
