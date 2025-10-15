function [xe,ye,se,xm,ym,sm]...
 ...
      = elm_line (N,ratio,x1,y1,x2,y2,sinit,Isym)

%-----------------------------------------
% FDLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement.
%----------------------------------------

%-----------------------------------------------------
% Disretization of a line segment into N elements
%
%  x1,y1: coordinates of the first point
%  x2,y2: coordinates of the last point
%  sinit: specified arc length at (x1, y1)
%
%  ratio: 
%
%  If Isym = 0, ratio of length of LAST to FIRST element
%  If Isym = 1, ratio of length of MID  to FIRST element
%  alpha: geometric factor ratio
%
%  se: polygonal arc length at the element end-nodes
%  sm: polygonal arc length at the element mid-nodes
%
%  xe,ye: end nodes
%  xm,ym: mid nodes
%-----------------------------------------------------

%------------
% one element
%------------

if(N==1) 

   xe(1) = x1;
   ye(1) = y1;

   xe(2) = x2; 
   ye(2) = y2;

   se(1) = sinit;
   se(2) = se(1)+sqrt( (x2-x1)^2+(y2-y1)^2);

%--------------------
% biased distribution
%--------------------

elseif(Isym==0)

      if(ratio==1.000)
        alpha  = 1.0;
        factor = 1.0/N;
      else
        texp   = 1.0/(N-1.0);
        alpha  = ratio^texp;
        factor = (1.0-alpha)/(1.0-alpha^N);
      end

      deltax = (x2-x1)*factor;  % x length of first element
      deltay = (y2-y1)*factor;  % y length of first element 

      xe(1) = x1;    % first point
      ye(1) = y1;
      se(1) = sinit;

      for i=2:N+1
        xe(i)  = xe(i-1)+deltax;
        ye(i)  = ye(i-1)+deltay;
        se(i)  = se(i-1)+sqrt(deltax^2+deltay^2);
        deltax = deltax*alpha;
        deltay = deltay*alpha;
      end

%----------------------
% symmetric distribution
%----------------------

elseif(Isym==1)

 xh = 0.5*(x1+x2);     % midpoint
 yh = 0.5*(y1+y2);

 if(N==2)

   xe(1) = x1;
   ye(1) = y1;

   xe(2) = xh;
   ye(2) = yh;

   xe(3) = x2;
   ye(3) = y2;

   se(1) = sinit;
   se(2) = se(1)+sqrt( (xh-x1)^2+(yh-y1)^2);
   se(3) = se(2)+sqrt( (x2-yh)^2+(y2-yh)^2);

 elseif(mod(N,2)==0) % even number of points

   Nh  = N/2;
   Nh1 = Nh+1;

      if(ratio==1.000)
        alpha  = 1.0;
        factor = 1.0/Nh;
      else
        texp   = 1.0/(Nh-1.0);
        alpha  = ratio^texp;
        factor = (1.0-alpha)/(1.0-alpha^Nh);
      end

      deltax = (xh-x1)*factor;  % x length of first element
      deltay = (yh-y1)*factor;  % y length of first element 

      xe(1) = x1;   % first point
      ye(1) = y1;
      se(1) = sinit;

      for i=2:Nh1
        xe(i)  = xe(i-1)+deltax;
        ye(i)  = ye(i-1)+deltay;
        se(i)  = se(i-1)+sqrt(deltax^2+deltay^2);
        deltax = deltax*alpha;
        deltay = deltay*alpha;
      end

      deltax = deltax/alpha;
      deltay = deltay/alpha;

      for i=Nh1+1:N+1
        xe(i)  = xe(i-1)+deltax;
        ye(i)  = ye(i-1)+deltay;
        se(i)  = se(i-1)+sqrt(deltax^2+deltay^2);
        deltax = deltax/alpha;
        deltay = deltay/alpha;
      end

   else

% odd number of points

     if(ratio==1.000)
        alpha  = 1.0;
        factor = 1.0/(N+1);
      else
        texp   = 2.0/(N-1.0);
        alpha  = ratio^texp;
        tmp1   = 0.50*(N+1.0);
        tmp2   = 0.50*(N-1.0);
        factor = (1.0-alpha)/(2.0-alpha^tmp1-alpha^tmp2);
     end

     deltax = (x2-x1)*factor;   % x length of first element
     deltay = (y2-y1)*factor;  % y length of first element 

     xe(1) = x1;    % first point
     ye(1) = y1;
     se(1) = sinit

     for i=2:(N+3)/2
        xe(i)  = xe(i-1)+deltax;
        ye(i)  = ye(i-1)+deltay;
        se(i)  = se(i-1)+sqrt(deltax^2+deltay^2);
        deltax = deltax*alpha;
        deltay = deltay*alpha;
     end

     deltax = deltax/(alpha^2);
     deltay = deltay/(alpha^2);

     for i=(N+5)/2:N+1
        xe(i)  = xe(i-1)+deltax;
        ye(i)  = ye(i-1)+deltay;
        se(i)  = se(i-1)+sqrt(deltax^2+deltay^2);
        deltax = deltax/alpha;
        deltay = deltay/alpha;
     end

end

%---
end
%---

%-------------------------------
% compute the element mid-points
%-------------------------------

  for i=1:N
    xm(i) = 0.5*(xe(i)+xe(i+1));
    ym(i) = 0.5*(ye(i)+ye(i+1));
    sm(i) = 0.5*(se(i)+se(i+1));
  end

%-----
% done
%-----

return
