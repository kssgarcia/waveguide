      subroutine prtcl_ax_vel
     +
     +    (Iflow
     +    ,X00,Y00
     +    ,Ux,Uy
     +    )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the lisencing agreement.
c----------------------------------------

c------------------------------------
c Evaluate the velocity at the point:
c             X00, Y00
c
c using the single-layer representation
c
c Capacity:
c --------
c
c   25 particles
c   64 elements along each particle
c
c--------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    NE(25),  Itp(25)
      Dimension xcntr(25),ycntr(25)
      Dimension axis1(25),axis2(25),tilt(25)

      Dimension xw(25,65),yw(25,65),tw(25,65)
      Dimension fx(25,64),fy(25,64)

      Dimension ux0(3200),uy0(3200)

c--------------
c common blocks
c--------------

      common/INTGR1/Nprtcl,NGL,NE,Itp,Ncl

      common/points/xw,yw,tw
      common/particles/xcntr,ycntr,axis1,axis2,tilt

      common/colloc1/fx,fy,ux0,uy0

      common/REAL1/visc,Uinf,wall,pg,RL,sc

      common/piii/pi,pi2,pi4,pi8

c--------
c prepare
c--------

      Ising = 0   ! will not desingularize the integrals

c-----------
c initialize
c-----------

      Ux = 0.0D0
      Uy = 0.0D0

      unused = 0.0D0

      j = 0         ! element counter

c---------------------------------
c Boundary integral representation
c---------------------------------

      Do K=1,Nprtcl

        Do L=1,NE(K)

        X1 = XW(K,L)
        Y1 = YW(K,L)
        T1 = TW(K,L)

        X2 = XW(K,L+1)
        Y2 = YW(K,L+1)
        T2 = TW(K,L+1)

        j = j+1

         call prtcl_ax_slp
     +
     +      (Iflow
     +      ,X00,Y00,unused
     +      ,X1,Y1,T1
     +      ,X2,Y2,T2
     +      ,NGL
     +      ,Ising
     +      ,Itp(k)
     +      ,xcntr(k),ycntr(k)
     +      ,axis1(k),axis2(k),tilt(k)
     +      ,Qxx,Qxy
     +      ,Qyx,Qyy
     +      )

        Ux = Ux-fx(K,L)*Qxx-fy(K,L)*Qxy
        Uy = Uy-fx(K,L)*Qyx-fy(K,L)*Qyy

        End Do

      End Do

c-------
c finish
c------

      Ux = Ux/pi8
      Uy = Uy/pi8

c----------------------
c add the incident flow
c----------------------

      If(Iflow.eq.1) then
         Ux = Ux+Uinf
      Else If(Iflow.eq.5) then
         Ux = Ux+0.25*pg/visc*(sc**2-Y00**2)
      End If

c-----
c Done
c-----

      Return
      End
