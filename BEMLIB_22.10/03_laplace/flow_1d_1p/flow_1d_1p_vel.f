      subroutine velocity 
     +
     +  (Xvel,Yvel
     +  ,shrt
     +  ,Uslip
     +  ,vel
     +  )

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-------------------------------------------------
c Computes the velocity at the point (Xvel, Yvel)
c using the bounday-integral representation
c-------------------------------------------------

      Implicit Double Precision (A-H,O-Z)

      Dimension NE(10),Itp(10)
      Dimension actis(10),xcntr(10),ycntr(10)

      Dimension xw(10,129),yw(10,129),tw(10,129)
      Dimension          dudn(10,128)

c---
c common blocks
c---

      common/VEL00/RL,visc,Iflow,NSG,NGL,NE,Itp
      common/VEL02/xw,yw,tw
      common/VEL03/actis,xcntr,ycntr
      common/VEL06/dudn

c-----
c constants
c-----

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi

c---
c initialize
c---

      Ising = 0    ! all elements are non-singular

      vel = 0.0D0

      unused = 0.0D0

c---------------------------------
c Boundary integral representation
c---------------------------------

c     test = 0.0D0

      j = 0         ! counter

      Do k=1,NSG

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

         call flow_1d_1p_slp
     +
     +      (Iflow
     +      ,RL
     +      ,Xvel,Yvel,unused
     +      ,X1,Y1,T1
     +      ,X2,Y2,T2
     +      ,NGL
     +      ,Ising
     +      ,Itp(k)
     +      ,rad,xcnt,ycnt
     +      ,QQQ
     +      )
c       write (6,*) dudn(k,l),QQQ

        vel = vel + dudn(k,l)*QQQ

c---
c       Dl  = Dsqrt((X2-X1)**2+(Y2-Y1)**2)
c       write (6,*) dudn(k,l),Dl,X1,X2,Y1,Y2
c       test = test + dudn(k,l)*Dl
c---

        End Do

      End Do

c     write (6,*) test
c     write (6,*) vel

c---
c complete the representation
c---

      vel = -vel/visc + shrt*Yvel + Uslip

c-----
c Done
c-----

      Return
      End
