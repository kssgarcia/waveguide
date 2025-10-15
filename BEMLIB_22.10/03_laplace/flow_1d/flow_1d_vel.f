      subroutine velocity 
     +
     +  (Xvel,Yvel
     +  ,dudn
     +  ,vel
     +  )

c=======================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=======================================

c-------------------------------------------------
c Compute the velocity at the point (Xvel, Yvel)
c using the bounday-integral representation
c-------------------------------------------------

      Implicit Double Precision (A-H,O-Z)

      Dimension NE(10),Itp(10)
      Dimension actis(10),xcntr(10),ycntr(10)

      Dimension xw(10,300),yw(10,300),tw(10,300)
      Dimension          dudn(10,300)

      Dimension u0(900)

c---
c common blocks
c---

      common/VEL00/Iflow,NSG,NGL,NE,Itp
      common/VEL02/xw,yw,tw,u0
      common/VEL03/actis,xcntr,ycntr

c----------
c constants
c----------

      pi  = 3.14159 265358
      pi2 = 2.0*pi

c-----------
c initialize
c-----------

      Ising = 0    ! all elements are non-singular

      vel = 0.0D0

      unused = 0.0

c---------------------------------
c Boundary integral representation
c---------------------------------

      j = 0         ! counter

      Do k=1,NSG

       rad  = actis(k)
       xcnt = xcntr(k)
       ycnt = ycntr(k)

        Do l=1,NE(k)

        X1 = XW(k,l)
        Y1 = YW(k,l)
        T1 = TW(k,l)

        X2 = XW(k,l+1)
        Y2 = YW(k,l+1)
        T2 = TW(k,l+1)

        j = j+1

        call flow_1d_sdlp
     +
     +      (Iflow
     +      ,Xvel,Yvel,unused
     +      ,X1,Y1,T1
     +      ,X2,Y2,T2
     +      ,NGL
     +      ,Ising
     +      ,Itp(k)
     +      ,rad,xcnt,ycnt
     +      ,QQQ
     +      ,WWW
     +      )

        vel = vel - dudn(k,l)*QQQ + u0(j)*WWW

        End Do

      End Do

c-----
c Done
c-----

      Return
      End
