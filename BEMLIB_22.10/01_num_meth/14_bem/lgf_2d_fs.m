function [G,Gx,Gy] = lgf_2d_fs (x,y,x0,y0,Iopt)

%-------------------------------------------
% Free-space Green's function:
%
%  G = -(1/2*pi) * lnr
%
% Iopt = 1: compute only the Green's function
%     ne 1: compute the Green's function
%           and the gradient
%
%-------------------------------------------

%-----------------
% Green's function
%-----------------

dx = x-x0;
dy = y-y0;
rs = dx*dx+dy*dy;
G  = - 0.5*log(rs)/(2*pi);

%--------------------------
% Green's function gradient
%--------------------------

if(Iopt~=1)
  den = rs*2*pi;
  Gx = - dx/den;
  Gy = - dy/den;
end

%-----
% done
%-----

return
