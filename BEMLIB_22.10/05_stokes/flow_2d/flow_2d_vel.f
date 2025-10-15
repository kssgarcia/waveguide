      subroutine velocity 
     +
     +   (X00,Y00
     +   ,Ux,Uy
     +   )

c-----------------------------------------
c FDLIB - BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c----------------------------------------
c  Computation of the velocity at the point
c  X00, Y00
c  using the boundary integral representation
c----------------------------------------

      Implicit Double Precision (A-H,O-Z)

      Dimension NE(10),Itp(10)
      Dimension actis(10),xcntr(10),ycntr(10)

      Dimension xw(10,200),yw(10,200),tw(10,200)

      Dimension fx(10,200),fy(10,200)
      Dimension ux0(900),uy0(900)

c---
c common blocks
c---

      common/xxx01/Iflow,NSG,NGL,NE,Itp
      common/xxx02/xw,yw,tw,fx,fy,ux0,uy0
      common/xxx03/actis,xcntr,ycntr
      common/xxx04/visc,shrt,delta

      common/flow_91/RL,Uslip

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi4 = 4.0D0*pi

c-----------
c initialize
c-----------

      Ising = 0    ! integrals are regular

      cf_slp = -1.0D0/(pi4*visc)
      cf_dlp =  1.0D0/ pi4

      Ux = 0.0D0
      Uy = 0.0D0

      unused = 0.0D0

      j = 0        ! element counter

c---------------------------------
c Boundary integral representation
c for the disturbance velocity
c---------------------------------

      Do K=1,NSG          ! loop over segments

       rad  = actis(k)
       xcnt = xcntr(k)
       ycnt = ycntr(k)

        Do L=1,NE(K)      ! loop over elements

        X1 = XW(K,L)
        Y1 = YW(K,L)
        T1 = TW(K,L)

        X2 = XW(K,L+1)
        Y2 = YW(K,L+1)
        T2 = TW(K,L+1)

        j = j+1

         call flow_2d_sdlp
     +
     +     (Iflow
     +     ,X00,Y00,unused
     +     ,X1,Y1,T1
     +     ,X2,Y2,T2
     +     ,NGL
     +     ,Ising
     +     ,Itp(k)
     +     ,rad,xcnt,ycnt
     +     ,Qxx,Qxy
     +     ,Qyx,Qyy
     +     ,Wxx,Wyx
     +     ,Wxy,Wyy
     +     )

        Ux = Ux + cf_slp*(fx(K,L)*Qxx+fy(K,L)*Qyx)
     +          + cf_dlp*(ux0(j) *Wxx+uy0(j) *Wyx)

        Uy = Uy + cf_slp*(fx(K,L)*Qxy+fy(K,L)*Qyy)
     +          + cf_dlp*(ux0(j) *Wxy+uy0(j) *Wyy)

        End Do

      End Do

c--------------------------
c add the incident velocity
c or the drift velocity
c--------------------------

      if(Iflow.eq.91.or.Iflow.eq.92) then

        Ux = Ux + shrt*Y00 + Uslip

      else

        Ux = Ux + shrt*Y00 + 0.5D0*delta/visc * Y00*Y00

      end if

c-----
c Done
c-----

      return
      end
