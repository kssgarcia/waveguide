function [f,Divcl,Divcint] ...
...
   = color(ax,bx,ay,by,Nx,Ny,NSG,shape)

%====================================
% FDLIB
%
% compute a mollified color function
% in a two-dimensional domain
%====================================

%---
% prepare
%---

L = bx-ax;
Dx = L/Nx;
Dy = (by-ay)/Ny;

Dx2 = 2.0*Dx;
Dy2 = 2.0*Dy;
Dxs = Dx*Dx;

%---
% parameters
%---

itermax = 128*128;
tol = 0.000000001;
relax = 1.0;

%--------
% interfacial markers on a circle
%--------

if(shape==1)

 radius = 0.10;
 centerx = 0.4;
 centery = 0.6;
 Dtheta = 2*pi/NSG;

 for k=1:NSG
  theta = (k-1)*Dtheta;
  XI(k) = centerx+radius*cos(theta);
  YI(k) = centery+radius*sin(theta);
  Dl(k) = radius*Dtheta;
  vnx(k) = cos(theta);  % normal vector
  vny(k) = sin(theta);
 end

 ftop = 0.0;
 fbot = 0.0;

%--------
% interfacial markers on a sinusoidal line
%--------

elseif(shape==2)

 Dxm = L/NSG;
 amp = 0.2;

 for k=1:NSG
  xaux = (k-0.5)*Dxm;
  XI(k) = xaux;
  YI(k) = 0.5+amp*sin(2*pi*xaux/L);
  vnx(k) = -amp*(2*pi/L)*cos(2*pi*xaux/L);
  vnm = sqrt(vnx(k)^2+1);
  vnx(k) = vnx(k)/vnm;
  vny(k) = 1.0/vnm;
  Dl(k) = Dxm*vnm;
 end

 ftop = 0.0;
 fbot = 1.0;

end

%--------
% project the gradient onto the grid
%--------

for i=1:Nx+2
 for j=1:Ny+1

   X(i,j) = (i-1)*Dx;
   Y(i,j) = (j-1)*Dy;
   cintx(i,j) = 0.0;
   cinty(i,j) = 0.0;
   Divcl(i,j) = 0.0;

%---
    for k=1:NSG     % loop over interfacial markers
%---

     yh = Y(i,j)-YI(k);

%---
      for ip=-1:1   % loop over periodic images
%---

       xh = X(i,j)-XI(k)+ip*L;
       fcx  = 0.0;
       fcxx = 0.0;
       fcyy = 0.0;
       if(abs(xh) <= 1.0*Dx2)
         fcx = (1.0+cos(pi*xh/Dx2))/(2.0*Dx2);
         fcxx = -pi*sin(pi*xh/Dx2)/(2.0*Dx2*Dx2);
       end
       fcy = 0.0;
       if(abs(yh) <= 1.0*Dy2) 
        fcy = (1.0D0+cos(pi*yh/Dy2))/(2.0*Dy2);
        fcyy = -pi*sin(pi*yh/Dy2)/(2.0*Dx2*Dy2);
       end

       fc = fcx*fcy*Dl(k);

       cintx(i,j) = cintx(i,j) + fc*vnx(k);
       cinty(i,j) = cinty(i,j) + fc*vny(k);
       Divcl(i,j) = Divcl(i,j) + Dl(k)*(fcxx*fcy*vnx(k)+ fcx*fcyy*vny(k));

       end
     end
%---
     end
 end

%-----------------------
% compute the divergence
%-----------------------

 Divcint(1:Nx+2,1:Ny+1) = 0.0;

 for i=2:Nx+1
  for j=2:Ny
   Divcint(i,j) = (cintx(i+1,j)-cintx(i-1,j))/Dx2 ...
                + (cinty(i,j+1)-cinty(i,j-1))/Dy2;
   end
 end

 % wrap

 for j=2:Ny
   Divcint(1,j) = Divcint(Nx+1,j);
   Divcint(Nx+2,j) = Divcint(2,j);
 end

%------
% solve
%------

[f,iter,Iflag] = pois_gs_dpr ...
...
   (Nx,Ny,Dx,Dy,Divcint,itermax,tol,relax,fbot,ftop);
%   (Nx,Ny,Dx,Dy,Divcl,itermax,tol,relax,fbot,ftop);

if(Iflag==0)
 disp "Solution did not coverge"
end

%---
% done
%---

return
