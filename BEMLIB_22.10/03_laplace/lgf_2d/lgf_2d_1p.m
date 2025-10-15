function [G,Gx,Gy] =  lgf_2d_1p ...
 ...
        (Iopt ...
        ,x,y ...
        ,x0,y0 ...
        ,RL ...
        )

%=======================================
% FDLIB, BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved
%
% This program is to be used only under the
% stipulations of the licensing agreement
%=======================================

%=======================================
% Singly periodic Green's function
% of Laplace's equation
% in two dimensions
%
% SYMBOLS:
% -------
%
% x, y :	coordinates of the field point
% x0,y0:	coordinates of the singular point
%
% RL:		period
%
% Iopt  = 1 produces G
% Iopt ne 1 produces G and its gradient (Gx,Gy)
%
% Notes:
% ----
%
% As y-y0 tends to infinity,
%
% G --> -(y-y0)/(2*RL)
%=======================================

%--------
% prepare
%--------

      wn = 2.0*pi/RL; %   wave number

      A = wn*(y-y0);
      B = wn*(x-x0);
      C = cosh(A)-cos(B);

%-----------------
% Green's function
%----------------- 

      G = -log(2.0*C)/(4*pi);

%--------------------------
% Green's function gradient
%--------------------------

      Gx = 0;
      Gy = 0;

      if(Iopt == 2)

      cf = -0.5/(RL*C);

      Gx = cf * sin(B);
      Gy = cf * sinh(A);

      end

%-----
% Done
%-----

      return
