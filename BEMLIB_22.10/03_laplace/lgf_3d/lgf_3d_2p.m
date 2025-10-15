function [G,Gx,Gy,Gz] = lgf_3d_2p ...
 ...
       (Iopt ...
       ,method ...
       ,x,y,z ...
       ,x0,y0,z0 ...
       ,a11,a12,a21,a22 ...
       ,b11,b12,b21,b22 ...
       ,ew,area ...
       ,Max1,Max2,Max3 ...
       )

%================================================
% FDLIB, BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved
%
% This program is to be used only under the
% stipulations of the licensing agreement
%================================================

%--------------------------------------
% Doubly periodic Green's function
% of Laplace's equation in three dimensions
%
% The point sources are located in a plane
% that is perpendicular to the z axis
%
% One of the point sources is located at (x0, y0, z0)
%
% The Green's function is computed by two methods
% according to the value of the switch "method"
%
% SYMBOLS:
% --------
%
% method = 1: Fourier series method
% method = 2: fast summation method of Hautman & Klein
%
% x, y, z:      coordinates of the field point
% x0, y0, z0:   coordinates of any one of the point sources
%
% G: Green's function
%
% Gx, Gy, Gz:  components of grad(G)
%
% Iopt =  1  compute only G
%      ne 1  compute G and grad(G)
%--------------------------------------

%----------
% constants
%----------

  pi2  = 2.0D0*pi;
  pi4  = 4.0D0*pi;
  srpi = sqrt(pi);

  toe =  3.0D0/8.0D0;
  tot =  3.0D0/2.0D0;
  foe = 15.0D0/8.0D0;

%--------
% prepare
%--------

  dxx = x-x0;
  dyy = y-y0;
  dzz = z-z0;

%-----------
% initialize
%-----------

  G  = 0.0D0; 
  Gx = 0.0D0;
  Gy = 0.0D0;
  Gz = 0.0D0;

%-----------------
% expedited method
%-----------------

if(method == 2)

%-----------------------------------------
% sum in real space of the 3D regularized
% Green's function over the 2D lattice
%-----------------------------------------

%-----
%  compute the accelerated sum of 1/r
%-----

  fc = 2.0/srpi;

  dz  = dzz;
  dzs = dz*dz;
  dzc = dzs*dz;
  dzq = dzs*dzs;

  for i1=-Max3:Max3
    for i2=-Max3:Max3

     dx = dxx - i1*a11 - i2*a21;
     dy = dyy - i1*a12 - i2*a22;

     %---
     if( (abs(dx)+abs(dy)+abs(dz)) > 0.00000001)

       sns = dx*dx+dy*dy;
       sn  = sqrt(sns);

       rns = sns+dzs;
       rn  = sqrt(rns);

       en  = dzs/sns;
       ens = en*en;

       G = G + 1.0/rn - (1.0-0.5*en+toe*ens)/sn;

      %  grad(G)

        if(Iopt>1) 
          rnc = rn*rns;
          snc = sn*sns;
          tmp = (1.0D0-1.5D0*en+foe*ens)/snc;
          Gx  = Gx-dx/rnc+dx*tmp;
          Gy  = Gy-dy/rnc+dy*tmp;
          Gz  = Gz-dz/rnc-(-1.0+tot*en)*dz/snc;
        end

     end
     %---

    end % for
  end % for

%-------------------------------------
% compute sums of 1/s**m in real space
%-------------------------------------

  sum0 = 0.0D0;
  sum1 = 0.0D0;
  sum2 = 0.0D0;

  if(Iopt >1 )
    sum0x = 0.0D0;
    sum0y = 0.0D0;
    sum1x = 0.0D0;
    sum1y = 0.0D0;
    sum2x = 0.0D0;
    sum2y = 0.0D0;
  end

  for i1=-Max1:Max1
    for i2=-Max1:Max1

        dx = dxx - i1*a11 - i2*a21;
        dy = dyy - i1*a12 - i2*a22;

        if((abs(dx)+abs(dy))>0.00000001) 

        sns = dx*dx+dy*dy;    % square of s_n
        sn  = sqrt(sns);      % s_n
        snc = sns*sn;         % thrid power of s_n
        snp = sns*snc;        % fifth power of s_n

        w  = ew*sn;
        ws = w*w;
        wq = ws*ws;
        we = wq*ws;

%---
% compute the error function
% erf(w)
%---

        T = 1.0/(1.0+0.5*w);
        erfcc =    ...         % complementary error function
        T*exp(-w*w-1.26551223+T*(1.00002368+T*(.37409196+ ...
        T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+...
        T*(1.48851587+T*(-.82215223+T*.17087277)))))))));

        erfc = 1.0D0-erfcc;

        expp = exp(-ws);

        h0 = erfc;
        h1 = h0- fc * w * expp*(1.0+2.0*ws);
        h2 = h0- fc * w * expp*(9.0+6.0*ws-4.0*wq+8.0*we)/9.0D0;

        sum0 = sum0+(1.0D0-h0)/sn;
        sum1 = sum1+(1.0D0-h1)/snc;
        sum2 = sum2+(1.0D0-h2)/snp;

%----------------
% compute grad(G)
%----------------

        if(Iopt>1) 

         snh = sns*snp;    % seventh power of s_n

         wq = ws*ws;
         dh0dw = fc*expp;
         dh1dw = dh0dw *  4.0*ws*(-1.0+ws);
         dh2dw = dh0dw * 16.0*wq*( 2.0-4.0*ws+wq)/9.0;

         tmp   = (1.0-h0+w*dh0dw)/snc;
         sum0x = sum0x - dx*tmp;
         sum0y = sum0y - dy*tmp;

         tmp   = (3.0-3.0*h1+w*dh1dw)/snp;
         sum1x = sum1x - dx*tmp;
         sum1y = sum1y - dy*tmp;

         tmp   = (5.0-5.0*h2+w*dh2dw)/snh;
         sum2x = sum2x - dx*tmp;
         sum2y = sum2y - dy*tmp;

        end
%------------------------- 
        end

       end % for
      end % for
%------------------------- 

      G = G + sum0 - 0.5D0*dzs*sum1  + toe*dzq*sum2;

%-------------------------        ! compute grad(G)
      if(Iopt>1)
        Gx = Gx + sum0x - 0.5*dzs*sum1x + toe*dzq*sum2x;
        Gy = Gy + sum0y - 0.5*dzs*sum1y + toe*dzq*sum2y;
        Gz = Gz - dz*sum1 + tot*dzc*sum2;
      end
%-------------------------

%-------------------------------
      end  % OF SUMMATION IN REAL SPACE
%-------------------------------

%-------------------------
% SUM IN WAVENUMBER SPACE
%-------------------------

      p  = 0.0D0;
      px = 0.0D0;
      py = 0.0D0;
      pz = 0.0D0;

      for i1=-Max2:Max2
       for i2=-Max2:Max2

        k1 = i1*b11 + i2*b21;
        k2 = i1*b12 + i2*b22;

        ks = k1*k1 + k2*k2;
        k  = sqrt(ks);

        if(k > 0.000001)  % skip the zero point

%-------------------------
      if(method == 1)  % straight Fourier series
%-------------------------

        rho = abs(dzz)*k;

        A = exp(-rho)/k;

        if(Iopt > 1)  % for grad(G)
         B =  k1*A;
         C =  k2*A;
         D = -k *A;
        end

%-------------------------
      elseif(method==2) % accelerated
%-------------------------

        w = k/(2.0D0*ew);
        T = 1.0D0/(1.0D0+0.5D0*w);          % complementary error function
        ERFCC= ...
        T*exp(-w*w-1.26551223+T*(1.00002368+T*(.37409196+ ...
        T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+ ...
        T*(1.48851587+T*(-.82215223+T*.17087277)))))))));

        fq0 =  erfcc/k;
        fq1 = -ks*fq0;
        fq2 = -ks*fq1/9.0D0;

        A = fq0 - 0.5D0 * dzs*fq1 + toe*dzq*fq2;

        if(Iopt > 1)     % for grad(G)
         B = k1*A;
         C = k2*A;
         D =  - dz*fq1 + tot*dzc*fq2;
        end

%-----------
      end
%-----------

      end

      arg = k1*dxx+k2*dyy;
      cs  = cos(arg);

      p = p + A*cs;

       if(Iopt > 1)     % for grad(G)
        sn = sin(arg);
        px = px - B*sn;
        py = py - C*sn;
        pz = pz + D*cs;
       end

   end  % for
   end  % for

%-------------------------------
% END OF SUM IN RECIPROCAL SPACE
%-------------------------------

%----------
% FINISH UP
%----------

%--------------------
% pure Fourier series
%-------------------

  if(method == 1)

        fc = 0.5/area;
        G = fc*p;
        G = G - fc*abs(dzz);     % regularize by adding a linear function

        if(Iopt>1)  % for grad(G)
         Gx = fc*px ;   % normalize
         Gy = fc*py;
         if(dzz<0) pz = -pz; end
         Gz = fc*pz - fc*dzz/abs(dzz); %  regularize 
        end

   end

%-----------------
% expedited method
%-----------------

 if(method == 2)

       prf = pi2/area;

       p = p - 1.0D0/(ew*srpi);  % contribution from zero wave number
       G = G + prf*p;             % add the sums in real and reciprocal space

       G = G/pi4;              % normalize

       if(Iopt>1) % for grad(G)

        Gx = Gx + prf*px;      % add the sums in real and reciprocal space
        Gy = Gy + prf*py;
        Gz = Gz + prf*pz;

        Gx = Gx/pi4;  % normalize
        Gy = Gy/pi4;
        Gz = Gz/pi4;

       end

 end

%---
% done
%---

return
