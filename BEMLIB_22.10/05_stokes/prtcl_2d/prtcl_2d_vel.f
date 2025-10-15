      subroutine velocity 
     +
     +    (Iflow
     +    ,X0,Y0
     +    ,Ux,Uy
     +    )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------
c Evaluate the velocity at the point:
c
c X0, Y0, 
c
c using the single-layer represenation
c
c Capacity:
c --------
c
c 72 particles
c 64 collocation points around each particle
c
c------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    NE(72),  Itp(72)
      Dimension axis1(72),axis2(72)
      Dimension xcntr(72),ycntr(72),tilt(72)

      Dimension xw(72,65),yw(72,65),tw(72,65)
      Dimension fx(72,64),fy(72,64)

      Dimension ux0(4608),uy0(4608)

c--------------
c common blocks
c--------------

      common/INTGR1/Nprtcl,NE,Itp,Ncl
      common/INTGR2/NGL

      common/REAL1/visc,shrt,delta,Uinf,wall

      common/CHANNELR/U1,U2,RL,h

      common/pax/axis1,axis2
      common/pap/xcntr,ycntr,tilt

      common/points/xw,yw,tw

      common/colloc1/ux0,uy0
      common/colloc6/fx,fy

      common/piii/pi,pi2,pi4

c--------
c prepare
c--------

      Ising = 0    ! will not desingularize the integrals

c-----------
c initialize
c-----------

      Ux = 0.0D0
      Uy = 0.0D0

      unused = 0.0D0

      Iopt = 1    ! will use the single-layer representation

      j = 0       ! element counter

c---------------------------------
c Boundary integral representation
c---------------------------------

      Do K=1,Nprtcl       ! loop over particles

        Do L=1,NE(K)      ! loop over elements

        X1 = XW(K,L)
        Y1 = YW(K,L)
        T1 = TW(K,L)

        X2 = XW(K,L+1)
        Y2 = YW(K,L+1)
        T2 = TW(K,L+1)

        j = j+1

         call prtcl_2d_sdlp
     +
     +      (Iopt
     +      ,Iflow
     +      ,X0,Y0,unused
     +      ,X1,Y1,T1
     +      ,X2,Y2,T2
     +      ,NGL
     +      ,Ising
     +      ,Itp(k)
     +      ,axis1(k),axis2(k)
     +      ,xcntr(k),ycntr(k),tilt(k)
     +      ,Qxx,Qxy
     +      ,Qyx,Qyy
     +      ,Wxx,Wyx
     +      ,Wxy,Wyy
     +      )

        Ux = Ux-fx(K,L)*Qxx-fy(K,L)*Qyx
c    +           +ux0(j)*Wxx+uy0(j)*Wyx
        Uy = Uy-fx(K,L)*Qxy-fy(K,L)*Qyy
c    +           +ux0(j)*Wxy+uy0(j)*Wyy

        End Do

      End Do

c-------

      Ux = Ux/(visc*pi4)
      Uy = Uy/(visc*pi4)

c----------------------
c Add the incident flow
c----------------------

c-----
      If(Iflow.eq.1.or.Iflow.eq.2          ! semi-infinite flow
     +             .or.Iflow.eq.3) then
c-----

       yref = Y0-wall
       Ux = Ux + shrt*yref + 0.5D0*delta*yref**2/visc

c-----
      Else If(Iflow.eq.3) then     ! channel flow
c-----

       Ux = Ux + U1+(U2-U1)*(Y0+h)/(2.0D0*h)
     +             +0.5D0*delta*Y0**2/visc

c-----
      Else If(Iflow.eq.10) then   ! periodic streaming flow
c-----

       Ux = Ux + Uinf

c-----------
      End If
c-----------

c-----
c Done
c-----

      Return
      End
