function [ber_0,bei_0,ber_0_p,bei_0_p] ...
...
      = ber_bei_0 (Iopt,X)

%===========================================
% FDLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licencing agreement
%===========================================

%-------------------------------------------
% Evaluation of the Kelvin functions:
%            ber_0, bei_0
%            and their derivatives
%
%  If Iopt.ne.1 derivatives are not computed
%-------------------------------------------

ber_0_p = 0;
bei_0_p = 0;

%---
 if(X<8.00)
%---

%---
% Use formulas (9.11.1) and (9.11.2) of A&S, p. 384
%---

  XH  = 0.5*X;
  Y   = X/8.0;
  Y2  = Y*Y;
  Y4  = Y2 *Y2;
  Y6  = Y4 *Y2;
  Y8  = Y6 *Y2;
  Y10 = Y8 *Y2;
  Y12 = Y10*Y2;
  Y14 = Y12*Y2;
  Y16 = Y14*Y2;
  Y18 = Y16*Y2;
  Y20 = Y18*Y2;
  Y22 = Y20*Y2;
  Y24 = Y22*Y2;
  Y26 = Y24*Y2;
  Y28 = Y26*Y2;

  ber_0 = 1.0 ...
          - 64.0        * Y4  + 113.77777774 * Y8 ...
          - 32.36345652 * Y12 +   2.64191397 * Y16 ...
          -  0.08349609 * Y20 +   0.00122552 * Y24 ...
          -  0.00000901 * Y28;
  bei_0 =    16.0        * Y2  - 113.77777774 * Y6 ...
          + 72.81777742 * Y10 -  10.56765779 * Y14 ...
          +  0.52185615 * Y18 -   0.01103667 * Y22 ...
          +  0.00011346 * Y26;
%---
   if(Iopt~=0)

%---
% use formulas (9.11.5) and (9.11.6) of A&S, p. 384
%---

    ber_0_p = X*(- 4.0        * Y2  + 14.22222222 * Y6 ...
                 - 6.06814810 * Y10 +  0.66047849 * Y14 ...
                 - 0.02609253 * Y18 +  0.00045957 * Y22 ...
                 - 0.00000394 * Y26 );
    bei_0_p = X*(    0.5              - 10.66666666 * Y4 ...
                  + 11.37777772 * Y8  -  2.31167514 * Y12 ...
                  +  0.14677204 * Y16 -  0.00379386 * Y20 ...
                  +  0.00004609 * Y24 );
   end

%---------
  else
%---------

  disp ' Sorry, ber_bei_0 not yet implemented for X>8'
  return

%---
 end
%---

return
