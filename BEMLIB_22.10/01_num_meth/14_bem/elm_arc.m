function [xe,ye,te,se,xm,ym,tm,sm]...
 ...
      = elm_arc (N,ratio,xcnt,ycnt ...
                ,radius,angle1,angle2 ...
                ,sinit,Isym)

%-----------------------------------------
% FDLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% this program is to be used only under the
% stipulations of the licensing agreement.
%----------------------------------------

%----------------------------------------------------------------
%  Disretization of a circular arc into N elements
%
%  SYMBOLS:
%  -------
%
%  xcnt, ycnt:	center of the arc
%  ratio:	radius of the arc
%  sinit: assigned arc length at the first point
%
%  ratio: if Isym = 0, ratio of length of LAST to FIRST element
%         if Isym = 1, ratio of length of MID  to FIRST element
%
%  xe, ye, te:	coordinates and polar angle for the end-nodes
%  xm, ym, tm:	coordinates and polar angle for the mid-nodes
%----------------------------------------------------------------

%------------
% one element
%------------

if(N==1)

 te(1) = angle1; 
 xe(1) = xcnt + radius*cos(te(1));
 ye(1) = ycnt + radius*sin(te(1));
 se(1) = sinit;

 te(2) = angle2;
 xe(2) = xcnt + radius*cos(te(2));
 ye(2) = ycnt + radius*sin(te(2));
 se(2) = se(1)+ radius*abs(te(2)-te(1));

%---------------------------
% non-symmetric distribution
%---------------------------

elseif(Isym==0)

   texp  = 1.0/(N-1.0);
   alpha = ratio^texp;
 
    if(abs(alpha-1.0)>0.0000001)
       factor = (1.0-alpha)/(1.0-alpha^N);
    else
       factor = 1.0/N;
    end

    deltat = (angle2-angle1)*factor;  % aperture of first element

      te(1) = angle1;                     % first point
      xe(1) = xcnt + radius*cos(te(1));
      ye(1) = ycnt + radius*sin(te(1));
      se(1) = sinit;

      for i=2:N+1
        te(i) = te(i-1)+deltat;
        xe(i) = xcnt + radius*cos(te(i));
        ye(i) = ycnt + radius*sin(te(i));
        se(i) = se(i-1)+radius*abs(deltat);
        deltat = deltat*alpha;
      end

%-----------------------
% symmetric distribution
%----------------------

  elseif(Isym==1)

  angleh = 0.5*(angle1+angle2);     % mid-point

%---
   if(N==2)
%---

   te(1) = angle1; 
   xe(1) = xcnt + radius*cos(te(1));
   ye(1) = ycnt + radius*sin(te(1));
   se(1) = sinit;

   te(2) = angleh;
   xe(2) = xcnt + radius*cos(te(2));
   ye(2) = ycnt + radius*sin(te(2));
   se(2) = se(1)+radius*abs(te(2)-te(1));

   te(3) = angle2;
   xe(3) = xcnt + radius*cos(te(3));
   ye(3) = ycnt + radius*sin(te(3));
   se(3) = se(2)+radius*abs(te(3)-te(2));

%---
  elseif(mod(N,2)==0) % even number of elements
%---

      Nh  = N/2;
      Nh1 = Nh+1;

      texp  = 1.0/(Nh-1.0);
      alpha = ratio^texp;

      if(abs(alpha-1.0)>0.0000001)
       factor = (1.0-alpha)/(1.0-alpha^Nh);
      else
       factor = 1.0/Nh;
      end

      deltat = (angleh-angle1)*factor;  % aperture of first element

      te(1) = angle1;                     % first point
      xe(1) = xcnt + radius*cos(te(1));
      ye(1) = ycnt + radius*sin(te(1));
      se(1) = sinit;

      for i=2:Nh1                          % up to mid-point
        te(i) = te(i-1)+deltat;
        xe(i) = xcnt + radius*cos(te(i));
        ye(i) = ycnt + radius*sin(te(i));
        se(i) = se(i-1)+radius*abs(deltat);
        deltat = deltat*alpha;
      end

      deltat = deltat/alpha;

      for i=Nh1+1:N+1                        % reflect
        te(i) = te(i-1)+deltat;
        xe(i) = xcnt + radius*cos(te(i));
        ye(i) = ycnt + radius*sin(te(i));
        se(i) = se(i-1)+radius*abs(deltat);
        deltat = deltat/alpha;
      end

%---
 else
%---

% odd number of elements

      texp  = 2.0/(N-1.0);
      alpha = ratio^texp;

      if(Dabs(alpha-1.0)>0.0000001)
       tmp1   = (N+1.0)/2.0;
       tmp2   = (N-1.0)/2.0;
       factor = (1.0-alpha)/(2.0-alpha^tmp1-alpha^tmp2);
      else
       factor = 1.0/N;
      end

      deltat = (angle2-angle1)*factor;  % aperture of first element

      te(1) = angle1;                     % first point
      xe(1) = xcnt + radius*cos(te(1));
      ye(1) = ycnt + radius*sin(te(1));
      se(1) = sinit;

      for i=2:(N+3)/2                      % up to mid point + 1
        te(i) = te(i-1)+deltat
        xe(i) = xcnt + radius*cos(te(i));
        ye(i) = ycnt + radius*sin(te(i));
        se(i) = se(i-1)+radius*abs(deltat);
        deltat = deltat*alpha;
      end

      deltat = deltat/(alpha^2);

      for i=(N+5)/2:N+1
        te(i) = te(i-1)+deltat;
        xe(i) = xcnt + radius*cos(te(i));
        ye(i) = ycnt + radius*sin(te(i));
        se(i) = se(i-1)+radius*abs(deltat);
        deltat = deltat/alpha;
      end

%---
   end
%---

%---
end
%---

%-----------------------
% compute the mid-points
%-----------------------

for i=1:N
  tm(i) = 0.5*(te(i)+te(i+1));
  sm(i) = 0.5*(se(i)+se(i+1));
  xm(i) = xcnt + radius*cos(tm(i));
  ym(i) = ycnt + radius*sin(tm(i));
end

%-----
% done
%-----

return
