function [Gxx,Gxy,Gxz ...
         ,Gyx,Gyy,Gyz ...
         ,Gzx,Gzy,Gzz ...
         ,px,py,pz ...
         ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz ...
         ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz ...
         ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz ...
         ] = ...
...
      sgf_3d_2p_w ...
...
       (Iopt ...
       ,x,y,z ...
       ,x0,y0,z0 ...
       ,wall ...
       ,a11,a12,a21,a22 ...
       ,Ns,Np ...
       )

%-----------------------------------------
% FDLIB, BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement.
%----------------------------------------

%-------------------------------------
% Periodic Green's function of 3D Stokes flow
% in a semi-infinite domain bounded by
% a plane wall located at x = wall
%
% The Green's function is computed by
% direct summation expedited by Aitken
% extrapolation
%-------------------------------------

%    if(Iopt==2)
%     Iopt = 1;
%     disp('Iopt reset to 1 in sgf_3d_2p_w')
%    end

%-----------
% initialize
%-----------

    Gxx = 0.0D0; Gxy = 0.0D0; Gxz = 0.0D0;
    Gyx = 0.0D0; Gyy = 0.0D0; Gyz = 0.0D0;
    Gzx = 0.0D0; Gzy = 0.0D0; Gzz = 0.0D0;

    px=0; py=0; pz=0;
    Txxx=0; Txxy=0; Txxz=0; Tyxy=0; Tyxz=0; Tzxz=0;
    Txyx=0; Txyy=0; Txyz=0; Tyyy=0; Tyyz=0; Tzyz=0;
    Txzx=0; Txzy=0; Txzz=0; Tyzy=0; Tyzz=0; Tzzz=0;

%---
% prepare
%---

      Ncut1 = Ns;
      Ncut2 = Ncut1*Np;
      Ncut3 = Ncut2*Np;

%---
% first run
%---

      savexx = 0.0D0;
      savexy = 0.0D0;
      savexz = 0.0D0;

      saveyx = 0.0D0;
      saveyy = 0.0D0;
      saveyz = 0.0D0;

      savezx = 0.0D0;
      savezy = 0.0D0;
      savezz = 0.0D0;

%---
% sum
%---

  for ip=-Ncut1:Ncut1
  for jp=-Ncut1:Ncut1

      x00 = x0;
      y00 = y0+ip*a11+jp*a21;
      z00 = z0+ip*a12+jp*a22;

      [Sxx,Sxy,Sxz ...
      ,Syx,Syy,Syz ...
      ,Szx,Szy,Szz ...
      ,px,py,pz ...
      ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz ...
      ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz ...
      ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz ...
      ] = ...
...
      sgf_3d_w ...
...
       (Iopt ...
       ,x,y,z ...
       ,x00,y00,z00 ...
       ,wall ...
       );

      Gxx = Gxx+Sxx;
      Gxy = Gxy+Sxy;
      Gxz = Gxz+Sxz;

      Gyx = Gyx+Syx;
      Gyy = Gyy+Syy;
      Gyz = Gyz+Syz;

      Gzx = Gzx+Szx;
      Gzy = Gzy+Szy;
      Gzz = Gzz+Szz;

      if((abs(ip)==Ncut1) | (abs(jp)==Ncut1)) 
       savexx = savexx+Sxx ;
       savexy = savexy+Sxy ;
       savexz = savexz+Sxz ;

       saveyx = saveyx+Syx ;
       saveyy = saveyy+Syy ;
       saveyz = saveyz+Syz ;

       savezx = savezx+Szx ;
       savezy = savezy+Szy ;
       savezz = savezz+Szz ;
      end

   end
   end

   gxxs1 = Gxx - 0.5*savexx ;
   gxys1 = Gxy - 0.5*savexy ;
   gxzs1 = Gxz - 0.5*savexz ;

   gyxs1 = Gyx - 0.5*saveyx ;
   gyys1 = Gyy - 0.5*saveyy ;
   gyzs1 = Gyz - 0.5*saveyz ;

   gzxs1 = Gzx - 0.5*savezx ;
   gzys1 = Gzy - 0.5*savezy ;
   gzzs1 = Gzz - 0.5*savezz ;

%---
% second run
%---

  savexx=0.0D0;
  savexy=0.0D0;
  savexz=0.0D0;

  saveyx=0.0D0;
  saveyy=0.0D0;
  saveyz=0.0D0;

  savezx=0.0D0;
  savezy=0.0D0;
  savezz=0.0D0;

  for ip=-Ncut2:Ncut2
  for jp=-Ncut2:Ncut2

      if((abs(ip)<=Ncut1) && (abs(jp)<=Ncut1)) 
%       Do nothing
      else

      x00 = x0 ;
      y00 = y0+ip*a11+jp*a21 ;
      z00 = z0+ip*a12+jp*a22 ;

      [Sxx,Sxy,Sxz ...
      ,Syx,Syy,Syz ...
      ,Szx,Szy,Szz ...
      ,px,py,pz ...
      ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz ...
      ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz ...
      ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz ...
      ] = ...
...
      sgf_3d_w ...
...
       (Iopt ...
       ,x,y,z ...
       ,x00,y00,z00 ...
       ,wall ...
       );

      Gxx = Gxx+Sxx;
      Gxy = Gxy+Sxy;
      Gxz = Gxz+Sxz;

      Gyx = Gyx+Syx;
      Gyy = Gyy+Syy;
      Gyz = Gyz+Syz;

      Gzx = Gzx+Szx;
      Gzy = Gzy+Szy;
      Gzz = Gzz+Szz;

      if((abs(ip)==Ncut2)|(abs(jp)==Ncut2))
       savexx=savexx+Sxx ;
       savexy=savexy+Sxy ;
       savexz=savexz+Sxz ;

       saveyx=saveyx+Syx ;
       saveyy=saveyy+Syy ;
       saveyz=saveyz+Syz ;

       savezx=savezx+Szx ;
       savezy=savezy+Szy ;
       savezz=savezz+Szz ;
      end

   end
  end
  end

  gxxs2 = Gxx - 0.5*savexx ;
  gxys2 = Gxy - 0.5*savexy ;
  gxzs2 = Gxz - 0.5*savexz ;

  gyxs2 = Gyx - 0.5*saveyx ;
  gyys2 = Gyy - 0.5*saveyy ;
  gyzs2 = Gyz - 0.5*saveyz ;

  gzxs2 = Gzx - 0.5*savezx ;
  gzys2 = Gzy - 0.5*savezy ;
  gzzs2 = Gzz - 0.5*savezz ;

%---
% third run
%---

  savexx=0.0D0;
  savexy=0.0D0;
  savexz=0.0D0;

  saveyx=0.0D0;
  saveyy=0.0D0;
  saveyz=0.0D0;

  savezx=0.0D0;
  savezy=0.0D0;
  savezz=0.0D0;

  for ip=-Ncut3:Ncut3
  for jp=-Ncut3:Ncut3

      if((abs(ip)<=Ncut2) && (abs(jp)<=Ncut2))
%       Do nothing
      else

      x00 = x0 ;
      y00 = y0+ip*a11+jp*a21 ;
      z00 = z0+ip*a12+jp*a22 ;

      [Sxx,Sxy,Sxz ...
      ,Syx,Syy,Syz ...
      ,Szx,Szy,Szz ...
      ,px,py,pz ...
      ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz ...
      ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz ...
      ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz ...
      ] = ...
...
      sgf_3d_w ...
...
       (Iopt ...
       ,x,y,z ...
       ,x00,y00,z00 ...
       ,wall ...
       );

      Gxx = Gxx+Sxx ;
      Gxy = Gxy+Sxy ;
      Gxz = Gxz+Sxz ;

      Gyx = Gyx+Syx ;
      Gyy = Gyy+Syy ;
      Gyz = Gyz+Syz ;

      Gzx = Gzx+Szx ;
      Gzy = Gzy+Szy ;
      Gzz = Gzz+Szz ;

      if((abs(ip)==Ncut3)|(abs(jp)==Ncut3))
       savexx=savexx+Sxx ;
       savexy=savexy+Sxy ;
       savexz=savexz+Sxz ;

       saveyx=saveyx+Syx ;
       saveyy=saveyy+Syy ;
       saveyz=saveyz+Syz ;

       savezx=savezx+Szx ;
       savezy=savezy+Szy ;
       savezz=savezz+Szz ;
      end

   end
  end
  end

      gxxs3 = Gxx - 0.5*savexx ;
      gxys3 = Gxy - 0.5*savexy ;
      gxzs3 = Gxz - 0.5*savexz ;

      gyxs3 = Gyx - 0.5*saveyx ;
      gyys3 = Gyy - 0.5*saveyy ;
      gyzs3 = Gyz - 0.5*saveyz ;

      gzxs3 = Gzx - 0.5*savezx ;
      gzys3 = Gzy - 0.5*savezy ;
      gzzs3 = Gzz - 0.5*savezz ;

%----
% extrapolate
%----

      Gxx = (gxxs3*gxxs1-gxxs2*gxxs2)/(gxxs3-2.0*gxxs2+gxxs1) ;
      Gxy = (gxys3*gxys1-gxys2*gxys2)/(gxys3-2.0*gxys2+gxys1) ;
      Gxz = (gxzs3*gxzs1-gxzs2*gxzs2)/(gxzs3-2.0*gxzs2+gxys1) ;

      Gyx = (gyxs3*gyxs1-gyxs2*gyxs2)/(gyxs3-2.0*gyxs2+gyxs1) ;
      Gyy = (gyys3*gyys1-gyys2*gyys2)/(gyys3-2.0*gyys2+gyys1) ;
      Gyz = (gyzs3*gyzs1-gyzs2*gyzs2)/(gyzs3-2.0*gyzs2+gyzs1) ;

      Gzx = (gzxs3*gzxs1-gzxs2*gzxs2)/(gzxs3-2.0*gzxs2+gzxs1) ;
      Gzy = (gzys3*gzys1-gzys2*gzys2)/(gzys3-2.0*gzys2+gzys1) ;
      Gzz = (gzzs3*gzzs1-gzzs2*gzzs2)/(gzzs3-2.0*gzzs2+gzzs1) ;

%----
% done
%----

  return
