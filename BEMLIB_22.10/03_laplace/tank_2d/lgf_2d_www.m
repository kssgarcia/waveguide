function [G,Gx,Gy]= lgf_2d_www ...
...
     (Iopt ...
     ,x,y ...
     ,x0,y0 ...
     ,wall1,wall2,wall3 ...
     )

%============================================
% FDLIB, BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved
%
% This program is to be used only under the
% stipulations of the licensing agreement
%============================================

%--------------------------------------------
% Neumman function of Laplace's equation in a 
% semi-infinite strip confined by three walls:
%
% a plane wall located at x = wall1
% a plane wall located at x = wall2
% a plane wall located at y = wall3
%
% Iopt =  1: compute only the Green's function 
%      ne 1: compute  the Green's function and gradient
%
% For formulae see Pozrikidis (1997, p. 364)
%--------------------------------------------

%--------
% prepare
%--------

 pi4 = 4.0*pi;

 h = wall2-wall1;
 h2 = 2.0*h;

 wn = 2*pi/h2;         % wave number

%--------------------------- 
% singularity and image system
% with respect to wall1
%--------------------------- 

 dx  = x-x0;
 dxi = x-2.0*wall1+x0;
 A   = wn*dx;
 Ai  = wn*dxi;
 dy = y-y0;
 B  = wn*dy;
 C  = cosh(B);
 D  = C-cos(A);
 Di = C-cos(Ai);

 G = -log(D*Di)/(4*pi);

%-----------------------------
% images with respect to wall3
%-----------------------------

 dy  = y-2.0*wall3+y0;
 BB  = wn*dy;
 CC  = cosh(BB);
 DD  = CC-cos(A);
 DDi = CC-cos(Ai);

 G = G - log(DD*DDi)/pi4;

%---------------------
% compute the gradient
%---------------------

  if(Iopt>1) 

    cf = -1.0D0/(2.0D0*h2);

    tmp  = 1.0D0/D +1.0D0/DD;
    tmpi = 1.0D0/Di+1.0D0/DDi;

    tlp  = 1.0D0/D +1.0D0/Di;
    tlpi = 1.0D0/DD+1.0D0/DDi;

    Gx = cf * ( sin(A)*tmp +  sin(Ai)*tmpi);
    Gy = cf * (sinh(B)*tlp + sinh(BB)*tlpi);

   end

%-----
% done
%-----

return
