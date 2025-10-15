function [al,be,ga] = abc ...
...

 (x1,y1,z1 ...
 ,x2,y2,z2 ...
 ,x3,y3,z3 ...
 ,x4,y4,z4 ...
 ,x5,y5,z5 ...
 ,x6,y6,z6 ...
 )

%========================================
% FDLIB, BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved
%
% This program is to be used only under
% the stipulations of the licensing agreement
%========================================

%--------------------------------------
% compute the parametric representation
% coefficients alpha, beta, gamma
%--------------------------------------

 d42 = sqrt( (x4-x2)^2 + (y4-y2)^2 + (z4-z2)^2 );
 d41 = sqrt( (x4-x1)^2 + (y4-y1)^2 + (z4-z1)^2 );
 d63 = sqrt( (x6-x3)^2 + (y6-y3)^2 + (z6-z3)^2 );
 d61 = sqrt( (x6-x1)^2 + (y6-y1)^2 + (z6-z1)^2 );
 d52 = sqrt( (x5-x2)^2 + (y5-y2)^2 + (z5-z2)^2 );
 d53 = sqrt( (x5-x3)^2 + (y5-y3)^2 + (z5-z3)^2 );

 al = 1.0/(1.0+d42/d41);
 be = 1.0/(1.0+d63/d61);
 ga = 1.0/(1.0+d52/d53);

%-----
% done
%-----

 return
