function [ew,area ...
         ,b11,b12,b21,b22] = lgf_3d_2p_ewald ...
     ...
       (a11,a12 ...
       ,a21,a22 ...
       )

%---------------------------------------------
% compute the reciprocal lattice base vectors
% and the optimal value of the parameter xi
%---------------------------------------------

%----------
% constants
%----------

  pi2  = 2.0*pi;
  srpi = sqrt(pi);

%---------------
% unit cell area
%---------------

  area = a11*a22-a21*a12;

%-----------------------------------------
% lattice base vectors in wave number space
%-----------------------------------------

  f    = pi2/area;
  b11  =  f*a22;
  b12  = -f*a21;
  b21  = -f*a12;
  b22  =  f*a11;

  ew = 0.5D0 * srpi/sqrt(area);

%-----
% done
%-----

return
