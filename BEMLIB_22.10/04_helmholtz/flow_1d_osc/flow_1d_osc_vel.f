      subroutine velocity 
     +
     +                (Xvel,Yvel
     +                ,delta
     +                ,velr,veli
     +                )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------
c Computes the real and imaginary parts
c of the velocity
c using the bounday-integral representation
c------------------------------------------

      Implicit Double Precision (A-H,O-Z)

      Dimension NE(10),Itp(10)
      Dimension actis(10),xcntr(10),ycntr(10)

      Dimension xw(10,200),yw(10,200),tw(10,200)
      Dimension dwdnr(10,200),dwdni(10,200)

      Dimension u0(500)

c--------------
c common blocks
c--------------

      common/VEL00/Iflow,NSG,NGL,NE,Itp
      common/VEL02/xw,yw,tw,u0
      common/VEL03/actis,xcntr,ycntr
      common/VEL05/dwdnr,dwdni

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi

c-----------
c initialize
c-----------

      Ising = 0    ! all elements are non-singular

      unused = 0.0D0

c---------------------------------
c Boundary integral representation
c---------------------------------

      velr = 0.0D0
      veli = 0.0D0

      j = 0         ! counter

      Do K=1,NSG

        If(Itp(k).eq.2) then
         rad  = actis(k)
         xcnt = xcntr(k)
         ycnt = ycntr(k)
        End If

        Do L=1,NE(K)

        X1 = XW(K,L)
        Y1 = YW(K,L)
        T1 = TW(K,L)

        X2 = XW(K,L+1)
        Y2 = YW(K,L+1)
        T2 = TW(K,L+1)

        j = j+1

         call flow_1d_osc_sdlp
     +
     +                   (Xvel,Yvel,unused
     +                   ,X1,Y1,T1
     +                   ,X2,Y2,T2
     +                   ,NGL
     +                   ,Ising
     +                   ,Itp(k)
     +                   ,delta
     +                   ,rad,xcnt,ycnt
     +                   ,QQQr,QQQi
     +                   ,WWWr,WWWi
     +                   )

        velr = velr - (dwdnr(k,l)*QQQr-dwdni(k,l)*QQQi)
     +              + u0(j)*WWWr

        veli = veli - (dwdnr(k,l)*QQQi+dwdni(k,l)*QQQr)
     +              + u0(j)*WWWi
        End Do

      End Do

c-----
c Done
c-----

      Return
      End
