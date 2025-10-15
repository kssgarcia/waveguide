      subroutine velocity 
     +
     +   (phi
     +   ,DphiDn0
     +   ,Iwall,wall
     +   ,X00,Y00
     +   ,Ux,Uy
     +   )

c==========================================
c FDLIB, BEMLIB, CFDLAB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c==========================================

c---------------------------------------
c Compute the velocity at a point using the 
c boundary-integral representation
c
c The velocity is computed by numerical differentiation
c of the potential, using centered differences.
c-----------------------------------------------------

      Implicit Double Precision (A-H,O-Z)

      Dimension    NE(10),Itp(10)
      Dimension actis(10),xcntr(10),ycntr(10)

      Dimension xw(10,129),yw(10,129),tw(10,129)
      Dimension phi(10,128)

      Dimension Dphidn0(1280)

      parameter (eps=0.01)

c--------------
c common blocks
c--------------

      common/xxx01/Iflow,NSG,NGL,NE,Itp
      common/xxx02/xw,yw,tw
      common/xxx03/actis,xcntr,ycntr
      common/xxx04/Vx,Vy,cr,xpv,ypv

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi

      eps2 = 2.0D0*eps   ! for centered differences

c------
c flags
c------

      Ising = 0

c-----------
c initialize
c-----------

      phi1 = 0.0D0
      phi2 = 0.0D0
      phi3 = 0.0D0
      phi4 = 0.0D0

      unused = 0.0

      X01 = X00-eps
      X02 = X00+eps

      Y01 = Y00-eps
      Y02 = Y00+eps

c---------------------------------
c Boundary integral representation
c---------------------------------

      j = 0         ! counter

      Do K=1,NSG

       rad  = actis(k)
       xcnt = xcntr(k)
       ycnt = ycntr(k)

        Do L=1,NE(K)

        X1 = XW(K,L)
        Y1 = YW(K,L)
        T1 = TW(K,L)

        X2 = XW(K,L+1)
        Y2 = YW(K,L+1)
        T2 = TW(K,L+1)

        j = j+1

         call body_2d_sdlp
     +
     +      (Iwall,wall
     +      ,X01,Y00,unused
     +      ,X1,Y1,T1
     +      ,X2,Y2,T2
     +      ,NGL
     +      ,Ising
     +      ,Itp(k)
     +      ,rad,xcnt,ycnt
     +      ,QQQ
     +      ,WWW
     +      )

        phi1 = phi1 - dphidn0(j)*QQQ + phi(k,l)*WWW

         call body_2d_sdlp
     +
     +      (Iwall,wall
     +      ,X02,Y00,unused
     +      ,X1,Y1,T1
     +      ,X2,Y2,T2
     +      ,NGL
     +      ,Ising
     +      ,Itp(k)
     +      ,rad,xcnt,ycnt
     +      ,QQQ
     +      ,WWW
     +      )

        phi2 = phi2 - dphidn0(j)*QQQ + phi(k,l)*WWW

         call body_2d_sdlp
     +
     +      (Iwall,wall
     +      ,X00,Y01,unused
     +      ,X1,Y1,T1
     +      ,X2,Y2,T2
     +      ,NGL
     +      ,Ising
     +      ,Itp(k)
     +      ,rad,xcnt,ycnt
     +      ,QQQ
     +      ,WWW
     +      )

        phi3 = phi3 - dphidn0(j)*QQQ + phi(k,l)*WWW

         call body_2d_sdlp
     +
     +      (Iwall,wall
     +      ,X00,Y02,unused
     +      ,X1,Y1,T1
     +      ,X2,Y2,T2
     +      ,NGL
     +      ,Ising
     +      ,Itp(k)
     +      ,rad,xcnt,ycnt
     +      ,QQQ
     +      ,WWW
     +      )

        phi4 = phi4 - dphidn0(j)*QQQ + phi(k,l)*WWW

        End Do

      End Do

c-------

      Ux = (phi2-phi1)/eps2 
      Uy = (phi4-phi3)/eps2

c--------------------------------
c add the incident streaming flow
c--------------------------------

      Ux = Ux + Vx 
      Uy = Uy + Vy

c------------------------------------------
c flow due to a point vortex of strength cr
c located at (Xpv,Ypv)
c------------------------------------------

      Xpvd = X00-Xpv
      Ypvd = Y00-Ypv
      Rpvds = Xpvd**2+Ypvd**2
      Ux = Ux - cr/pi2 * Ypvd/Rpvds 
      Uy = Uy + cr/pi2 * Xpvd/Rpvds

c------------------------------------------------
c image point vortex
c with respect to a flat wall located at y = wall
c------------------------------------------------

      If(Iwall.eq.1) then

        Xpvdi  = X00-Xpv
        Ypvdi  = Y00+Ypv-2.0*wall
        Rpvdsi = Xpvdi**2+Ypvdi**2
        Ux = Ux + cr/pi2 * Ypvdi/Rpvdsi 
        Uy = Uy - cr/pi2 * Xpvdi/Rpvdsi

      End If

c-----
c Done
c-----

      Return
      End
