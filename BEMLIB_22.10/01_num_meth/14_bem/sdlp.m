function [SLP, DLP] = sdlp ...
      ...
   (x0,y0,t0 ...    % evaluation point
   ,x1,y1,t1 ...    % first element point
   ,x2,y2,t2 ...    % second element point
   ,NGL ...         % Gauss--Legendre quadrature order
   ,Ising ...       % element singularity index
   ,Itype ...       % element type
   ,rad,xcnt,ycnt)  % arc radius and center

%======================================================
%
% FDLIB
%
% Compute the single-layer and double-layer potential 
% along a straight segment or a circular arc
%
% SYMBOLS:
% -------
%
% SLP: single-layer potential
% DLP: double-layer potential
%
% Ising = 0 : non-singular element
%         1 : singular element
%======================================================

%--------
% prepare
%--------

Iopt = 2;   % for the Green's function

SLP = 0.0;
DLP = 0.0;

[ZZ,WW] = gauss_leg(NGL);

%---------------------------
% prepare for the quadrature
%---------------------------

if(Itype==1) % straight segments

  xM = 0.5D0*(x2+x1);
  xD = 0.5D0*(x2-x1);
  yM = 0.5D0*(y2+y1);
  yD = 0.5D0*(y2-y1);
  DR = sqrt(xD*xD+yD*yD);
  vnx = -yD/DR;        % unit normal vector
  vny =  xD/DR;

else % circular arcs

  tM = 0.5*(t2+t1);
  tD = 0.5*(t2-t1);
  DR = rad*abs(tD);
  ornt = 1.0;         % orientation index
  if(tD<0) ornt = -1.0; end

end

%---
% loop over Gaussian points
%---

for i=1:NGL

   if(Itype==1)   % straight segments
     x = xM + xD*ZZ(i);
     y = yM + yD*ZZ(i);
    else          % circular arcs
     t  = tM + tD*ZZ(i);
     cs = cos(t);
     sn = sin(t);
     x  = xcnt + rad*cs;
     y  = ycnt + rad*sn;
     vnx = -cs*ornt;  % unit normal vector
     vny = -sn*ornt;  % when arc is clockwise,
    end               % normal vector points toward the center

    [G,dGdx,dGdy] = lgf_2d_fs (x,y,x0,y0,Iopt);

%---
% treat the slp singularity
%---

  if(Ising==1)
   if(Itype==1) Dists = (x-x0)^2+(y-y0)^2; end
   if(Itype==2) Dists = rad^2*(t0-t)^2; end
   G  = G + log(Dists)/(4*pi);
  end

  SLP = SLP + G*WW(i);
  DLP = DLP + (dGdx*vnx+dGdy*vny)*WW(i);

end

%---
% finish up
%---

SLP = SLP * DR;
DLP = DLP * DR;

%---
% add the slp singularity back to the slp
%---

if(Ising==1) 
  SLP = SLP - 2.0*DR*(log(DR)-1.0)/(2*pi);
end

%----------------------------------
% Analytical integration of the dlp
% for the free space GF
% over the singular elements
%
% Note that the dlp of the free-space
% Green's function vanishes over
% singular straight segments
%----------------------------------

if(Ising==1)
  if(Itype==1)             % straight segments
    DLP = 0.0;
  elseif(Itype==2)         % circular arcs
    DLP = (t2-t1)/(4*pi);  % independent of t0 !
  end
end

%-----
% done
%-----

return
