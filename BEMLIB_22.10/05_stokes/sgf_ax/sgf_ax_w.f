      subroutine sgf_ax_w
     +
     +  (Iopt
     +  ,X0,Y0
     +  ,X,Y
     +  ,wall
     +  ,SXX,SXY
     +  ,SYX,SYY
     +  ,QXXX,QXXY,QXYX,QXYY
     +  ,QYXX,QYXY,QYYX,QYYY
     +  ,PXX,PXY,PYX,PYY
     +  ,Iaxis
     +  )

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c---------------------------------------------------
c  Axisymmetric Green's function of Stokes flow
c  and kernel of the double-layer potential
c  for flow bounded by a plane wall located at x = wall
c
c  (x,   y) are the x and s coordinates of the singularity
c  (x0, Y0) are the x and s coordinates of the field point
c
c  Sij     Green's function
c  Qijk    Double-layer kernel
c
c  Iaxis = 0  field point is off-axis s = 0
c        = 1  field point is on-axis
c
c  wall:   x location of the wall
c
c  h:  distance of field point from wall
c
c  Use with the boundary integral equation:
c  ---------------------------------------
c
c  Let b be the strength of the point-force ring located at x;
c  then the induced velocity field is:
c
c  ux(x0) = Sxx(x,x0) * bx + Sxs(x,x0) * bs
c  us(x0) = Ssx(x,x0) * bx + Sss(x,x0) * bs
c
c  The kernel of the axisymmetric double-layer potential is:
c
c   Idlpx(x0) = ux * ( Qxxx * vnx + Qxxs * vns)
c             + us * ( Qxsx * vnx + Qxss * vns)
c
c   Idlps(x0) = ux * ( Qsxx * vnx + Qsxs * vns)
c             + us * ( Qssx * vnx + Qsss * vns)
c
c  The arguments of Qxxx are (x,x0)
c  This is the flow due to a ring distribution of stresslets
c
c  Pij is used for the desingularization of the dlp
c
c------------
c
c  Iopt = 1 produces only the Green's function
c  Iopt = 2 produces the Green's function and the stress tensor
c
c-------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Parameter (eps=0.00000001)

c----------
c constants
c----------

      pi  = 3.14159 265358 97932384 D0
      pi2 = 2.0D0*pi

c----------------------
c field point on axis ?
c----------------------

      Iaxis = 0
      if(Y0.lt.eps) Iaxis = 1

c---------------
c prepare to run
c---------------

      H    = X0-wall
      HH   = 2.0D0*H
      H2   = H*H
      H2H2 = H2+H2

      Y2  = Y*Y
      Y02 = Y0*Y0
      YY2 = Y2+Y02
      YYP = Y*Y0
      DY  = Y-Y0
      DY2 = DY*DY

      if(Iopt.ne.1) then
        Y3  = Y*Y2
        Y66 = 6.0D0*Y
        Y03 = Y0*Y02
      end if
 
      if(Iaxis.eq.0) then
        Y4  = Y2*Y2
        Y04 = Y02*Y02
        YYR = dsqrt(YYP)
        YY3 = YYR*YYR*YYR
        YY5 = YYR*YYR*YY3
        SY  = Y+Y0
        SY2 = SY*SY
      end if
 
      DX = X-X0

c-----------------------------
c will pass twice:
c
c  first for the primary ring
c  second for the image ring
c-----------------------------

c---
      Do Ipass=1,2
c---

      DX2  = DX*DX
      DX4  = DX2*DX2
      DR2  = DX2+DY2
      DR   = dsqrt(DR2)
      DXYY = DX2+YY2

      if(Iopt.ne.1) then
        Y6DX  = Y66*DX
        Y6DX2 = Y66*DX2
        DX3   = DX*DX2
        Y6DX3 = Y66*DX3
      end if
 
c----------------------
      if(Iaxis.eq.0) then  ! point off axis

        FC1  = 4.0D0*YYP
        FC2  = DX2+SY2
        RK2  = FC1/FC2
        RK   = dsqrt(RK2)
        RK3  = RK2*RK
        RK4  = RK2*RK2
        RK5  = RK4*RK
        RK2P = 1.0D0-RK2

        call ell_int (RK2,F,E)

        RJ10 = 2.0D0*RK*F/YYR
        RJ11 = RK*(DXYY*F-(DX2+SY2)*E)/YY3
        RL30 = E/RK2P
        FCTR = RK3/(2.0D0*YY3)
        RJ30 = FCTR*RL30
c       RJ30 = 2.0D0*RK*E/(YYR*DR2)
        RJ31 = RK*(-F+DXYY*E/DR2)/YY3
        RJ32 = RK*(-DXYY*F+(DX4+2.0*DX2*YY2+Y4+Y04)*E/DR2)/YY5
 
      if(Iopt.ne.1.or.Ipass.eq.2) then

        RK6  = RK4*RK2
c       RK8  = RK6*RK2
        RL50 =  (2.0D0*(2.0D0-RK2)*RL30-F)/(3.0D0*RK2P)
        RL52 =  (RL50 - RL30) / RK2
        RL54 =  (F-RL50+2.0D0*RK2*RL52) / RK4
        RL56 = -(E-RL50+3.0*RK2*RL52-3.0*RK4*RL54)/ RK6
c       PREP = ( 2.0*(2.0-RK2)*E - RK2P*RK2P*F ) / 3.0
c       RL58 = - (PREP - RL50 + 4.0*RK2*RL52 - 6.0*RK4*RL54
c    +                        + 4.0*RK6*RL56 ) / RK8
        FCTR = RK5/(8.0D0*YY5)
        RJ50 = FCTR * RL50
        RJ51 = FCTR * (2.0D0*RL52-RL50)
        RJ52 = FCTR * (4.0D0*(RL54-RL52)+RL50)
        RJ53 = FCTR * (8.0D0*RL56 - 12.0D0*RL54
     +                +6.0D0*RL52 - RL50)
c       RJ54 = FCTR * (16.0D0*RL58 - 32.0D0*RL56 + 24.0D0*RL54
c    +                 -8.0D0*RL52 + RL50)

      end if

c---
      else   ! point on axis
c---

        DR3  = DR2*DR
        DR5  = DR3*DR2
        RJ10 = pi2/DR
        RJ11 = 0.0D0
        RJ30 = PI2/DR3
        RJ31 = 0.0D0
        RJ32 = pi/DR3

        if(Iopt.ne.1.or.Ipass.eq.2) then
          RJ50 = pi2/DR5
          RJ51 = 0.0D0
          RJ52 = pi/DR5
          RJ53 = 0.0D0
C         RJ54 = 3.0D0*pi/(4.0*DR5)
        end if

c---
      end if 
c---

c-----------------
c Green's function
c-----------------
 
      SXX = Y *    (  RJ10+DX2*RJ30)
      SXY = Y * DX*(Y*RJ30-Y0 *RJ31)
      SYX = Y * DX*(Y*RJ31-Y0 *RJ30)
      SYY = Y *    (  RJ11+YY2*RJ31-YYP*(RJ30+RJ32))

      if(Iopt.ne.1) then

        QXXX = - Y6DX3 * RJ50
        QXXY = - Y6DX2 * (Y*RJ50 - Y0*RJ51)
        QXYX =   QXXY
        QXYY = - Y6DX  * (Y02*RJ52+Y2*RJ50 - 2.0*YYP*RJ51)
        QYXX = - Y6DX2 * (Y  *RJ51-Y0*RJ50)
        QYXY = - Y6DX  * ((Y02+Y2)*RJ51 - YYP*(RJ52+RJ50))
        QYYX = QYXY
        QYYY = - Y66 * ( Y3     *  RJ51 
     +                 - Y03    *  RJ52 
     +                 - Y*YYP  * (RJ50+2.0*RJ52)
     +                 + Y0*YYP * (RJ53+2.0*RJ51) )

        if(Ipass.eq.1) then
         PXX = QYXX
         PXY = QYXY
         PYX = - Y6DX * (Y2*RJ52 - 2.0*YYP*RJ51 + Y02*RJ50)
         PYY = - Y66  * (  Y0*YYP * (2.0*RJ52+RJ50)
     +                   - Y03    * RJ51
     +                   - Y *YYP * (2.0*RJ51+RJ53)
     +                   + Y3     * RJ52  )
         end if

      end if

c     write (6,*) SXX,SXY,SYX,SYY
c     write (6,*) QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY
c     stop
 
c--------------------
c IMAGE SINGULARITIES
c--------------------

      if(Ipass.eq.1) then
        SXXSAVE = SXX
        SXYSAVE = SXY
        SYXSAVE = SYX
        SYYSAVE = SYY
        if(Iopt.ne.1) then
          QXXXSAVE = QXXX
          QXXYSAVE = QXXY
          QXYXSAVE = QXYX
          QYXXSAVE = QYXX
          QXYYSAVE = QXYY
          QYXYSAVE = QYXY
          QYYXSAVE = QYYX
          QYYYSAVE = QYYY
        end if
        DX = X+X0-2.0D0*wall
      else
        SXX = SXXSAVE - SXX
        SXY = SXYSAVE - SXY
        SYX = SYXSAVE - SYX
        SYY = SYYSAVE - SYY
        if(Iopt.ne.1) then
          QXXX = QXXXSAVE - QXXX
          QXXY = QXXYSAVE - QXXY
          QXYX = QXYXSAVE - QXYX
          QYXX = QYXXSAVE - QYXX
          QXYY = QXYYSAVE - QXYY
          QYXY = QYXYSAVE - QYXY
          QYYX = QYYXSAVE - QYYX
          QYYY = QYYYSAVE - QYYY
        end if
      end if

c---
      end do   ! of Ipass
c---
 
c----------------------------------
c  SOURCE-DIPOLE AND STOKES-DOUBLET
c----------------------------------

      DDXX = Y *        (-RJ30+3.0*DX2*RJ50)
      DDXY = Y * 3.0*DX*( Y*RJ50-Y0*RJ51)
      DDYX = Y * 3.0*DX*(-Y*RJ51+Y0*RJ50)
      DDYY = Y *     (RJ31-3.0*YY2*RJ51+3.0*YYP*(RJ50+RJ52))
 
      EEXX = DX*DDXX
      EEXY = DX*DDXY + Y*(Y0*RJ31-Y*RJ30)
      EEYX = DX*DDYX + Y*(Y0*RJ30-Y*RJ31)
      EEYY = DX*DDYY
 
      SXX = SXX + H2H2*DDXX - HH*EEXX
      SXY = SXY + H2H2*DDXY - HH*EEXY
      SYX = SYX + H2H2*DDYX - HH*EEYX
      SYY = SYY + H2H2*DDYY - HH*EEYY

c---
      if(Iopt.ne.1) then
c---

      if(Iaxis.eq.0) then
        RL70 =  (0.80*(2.0-RK2)*RL50 - 0.60*RL30) / RK2P
        RL72 =  (RL70-RL50) / RK2
        RL74 =  (RL30-RL70+2.0D0*RK2*RL72) / RK4
        RL76 = -(F-RL70+3.0D0*RK2*RL72-3.0D0*RK4*RL74) / RK6
        RK7  = RK2*RK5
        YY7  = YYR*YYR*YY5
        FCT7 = RK7/(32.0D0*YY7)
        RJ70 = FCT7* RL70
        RJ71 = FCT7*(2.0D0*RL72-RL70)
        RJ72 = FCT7*(4.0D0*(RL74-RL72)+RL70)
        RJ73 = FCT7*(8.0D0*RL76-12.0D0*RL74+6.0D0*RL72-RL70)
      else
        DR7  = DR5*DR2
        RJ70 = pi2/DR7
        RJ71 = 0.0D0
        RJ72 = pi/DR7
        RJ73 = 0.0D0
      end if
c---
 
      DQXXX = - Y6DX * (-3.0*RJ50 + 5.0*DX2*RJ70)
      DQXXY = - Y66  * (-Y*RJ50+Y0*RJ51 - 5.0*DX2*(-Y*RJ70+Y0*RJ71))
      DQXYX =   DQXXY
      DQXYY = - Y6DX * (-RJ50 + 5.0*(Y2*RJ70-2.0*YYP*RJ71+Y02*RJ72))
      DQYXX =   Y66  * (-Y*RJ51+Y0*RJ50 - 5.0*DX2*(-Y*RJ71+Y0*RJ70))
      DQYXY =   Y6DX * (-RJ51 + 5.0*YY2*RJ71
     +                        - 5.0*YYP*(RJ70+RJ72))
      DQYYX = DQYXY
      DQYYY =    Y66*(-3.0D0*Y*RJ51+Y0*(RJ50+2.0D0*RJ52)
     +                +5.0D0*Y*Y02*RJ73
     +                +5.0D0*Y*(Y2+2.0*Y02)*RJ71
     +                -5.0D0*Y0*(2.0*Y2+Y02)*RJ72
     +                -5.0D0*Y2*Y0*RJ70)              
 
      FQXXX =   0.0D0
      FQXXY =   Y6DX * (Y*RJ50 - Y0*RJ51)
      FQXYX =   FQXXY
      FQXYY = - Y66  * (DX2*RJ50-Y2*RJ50-Y02*RJ52+2.0*YYP*RJ51)
      FQYXX =   Y6DX * (Y*RJ51 - Y0*RJ50)
      FQYXY =   0.0D0
      FQYYX =   0.0D0
      FQYYY =   Y6DX * (Y*RJ51 - Y0*RJ50)
 
      EQXXX = DX*DQXXX + FQXXX
      EQXXY = DX*DQXXY + FQXXY
      EQXYX = DX*DQXYX + FQXYX
      EQYXX = DX*DQYXX + FQYXX
      EQXYY = DX*DQXYY + FQXYY
      EQYXY = DX*DQYXY + FQYXY
      EQYYX = DX*DQYYX + FQYYX
      EQYYY = DX*DQYYY + FQYYY
 
      QXXX = QXXX + H2H2*DQXXX - HH*EQXXX
      QXXY = QXXY + H2H2*DQXXY - HH*EQXXY
      QXYX = QXYX + H2H2*DQXYX - HH*EQXYX
      QYXX = QYXX + H2H2*DQYXX - HH*EQYXX
      QXYY = QXYY + H2H2*DQXYY - HH*EQXYY
      QYXY = QYXY + H2H2*DQYXY - HH*EQYXY
      QYYX = QYYX + H2H2*DQYYX - HH*EQYYX
      QYYY = QYYY + H2H2*DQYYY - HH*EQYYY

c---
      end if
c---
 
c-----
c Done
c-----

      return
      end
