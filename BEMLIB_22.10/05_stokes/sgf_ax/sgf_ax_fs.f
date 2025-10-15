      subroutine sgf_ax_fs
     +
     + (Iopt
     + ,X0,Y0
     + ,X,Y
     + ,SXX,SXY
     + ,SYX,SYY
     + ,QXXX,QXXY,QXYX,QXYY
     + ,QYXX,QYXY,QYYX,QYYY
     + ,PXX,PXY,PYX,PYY
     + ,Iaxis
     + )

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c================================================
c  Axisymmetric Green's function of Stokes flow
c  in free space
c
c  Let b be the strength of the point-force ring located at x;
c  then the induced velocity field is:
c
c  ux(x0) = Sxx(x0,x) * bx + Sxs(x0,x) * bs
c  us(x0) = Ssx(x0,x) * bx + Sss(x0,x) * bs
c
c  The kernel of the axisymmetric double-layer potential is:
c
c   Idlpx(x0) = ux * ( Qxxx * vnx + Qxxs * vns)
c             + us * ( Qxsx * vnx + Qxss * vns)
c
c   Idlps(x0) = ux * ( Qsxx * vnx + Qsxs * vns)
c             + us * ( Qssx * vnx + Qsss * vns)
c
c  arguments of Qxxx are (x,x0)
c
c  This is the flow due to a ring distribution of stresslets
c
c  Pij is used for the desingularization of the dlp
c
c------------
c
c  Iopt = 1 produces only the Green's function
c  Iopt = 2 produces the Green's function and the stress tensor
c
c================================================

      Implicit Double Precision (a-h,o-z)

      Parameter (eps=0.00000001)

c----------
c constants
c----------

      pi  = 3.1415 92653 58979 32384 D0
      pi2 = 2.0D0*pi

c------------------------
c field point on x axis ?
c------------------------

      Iaxis = 0

      if(y0.lt.eps) Iaxis = 1

c-------
c prepare
c-------

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
 
      if(Iaxis.eq.0) then   ! off the axis
        Y4  = Y2*Y2
        Y04 = Y02*Y02
        YYR = sqrt(YYP)
        YY3 = YYR**3
        YY5 = YYR**5
        SY  = Y+Y0
        SY2 = SY*SY
      end if

      DX   = X-X0
      DX2  = DX*DX
      DX4  = DX2*DX2
      DR2  = DX2+DY2
      DR   = sqrt(DR2)
      DXYY = DX2+YY2

      if(Iopt.ne.1) then
        Y6DX  = Y66*DX
        Y6DX2 = Y66*DX2
        DX3   = DX*DX2
        Y6DX3 = Y66*DX3
      end if

c------------------------
      if(Iaxis.eq.0) then       !  point x0 off the axis
c------------------------

        FC1  = 4.0D0*YYP
        FC2  = DX2+SY2
        RK2  = FC1/FC2
        RK   = sqrt(RK2)
        RK3  = RK2*RK
        RK4  = RK2**2
        RK5  = RK4*RK
        RK2P = 1.0D0-RK2

        call ell_int (RK2,F,E)

        RI10 = 2.0D0*RK*F/YYR
        RI11 = RK*(DXYY*F-(DX2+SY2)*E)/YY3
        RI30 = 2.0D0*RK*E/(YYR*DR2)
        RI31 = RK*(-F+DXYY*E/DR2)/YY3
        RI32 = RK*(-DXYY*F+(DX4+2.0*DX2*YY2+Y4+Y04)*E/DR2)/YY5

        if(Iopt.ne.1) then

        RL10 = F
        RL30 = E/RK2P
        RK6  = RK4*RK2
c       RK8  = RK6*RK2
        RL50 =  (2.0D0 * (2.0D0 - RK2)*RL30 - F) / (3.0D0*RK2P)
        RL52 =  (RL50-RL30) / RK2
        RL54 =  (F-RL50+2.0*RK2*RL52) / RK4
        RL56 = -(E-RL50+3.0*RK2*RL52-3.0*RK4*RL54)/ RK6
c       PREP = ( 2.0*(2.0-RK2)*E - RK2P*RK2P*F ) / 3.0D0
c       RL58 = - (PREP - RL50 + 4.0D0 * RK2*RL52 - 6.0D0 * RK4*RL54
c    +                        + 4.0D0 * RK6*RL56 ) / RK8
        FCTR = RK5/(8.0*YY5)
        RI50 = FCTR * RL50
        RI51 = FCTR * (2.0D0*RL52-RL50)
        RI52 = FCTR * (4.0D0*(RL54-RL52)+RL50)
        RI53 = FCTR * (8.0D0*RL56 - 12.0D0*RL54
     +                +6.0D0*RL52 - RL50)
c       RI54 = FCTR * (16.0*RL58 - 32.0*RL56 + 24.0*RL54
c    +                 -8.0*RL52 + RL50)
        End If

c-----------
        else               !  point x0 on the axis
c-----------

        DR3  = DR2*DR
        DR5  = DR3*DR2
        RI10 = PI2/DR
        RI11 = 0.0D0
        RI30 = PI2/DR3
        RI31 = 0.0D0
        RI32 = PI/DR3

        if(Iopt.ne.1) then
         RI50 = PI2/DR5
         RI51 = 0.0D0
         RI52 = PI/DR5
         RI53 = 0.0D0
c        RI54 = 3.0*PI/(4.0*DR5)
        end if

c-------------
        end if
c-------------

c---
c build the Green's function
c---

      SXX = Y *    (  RI10 + DX2*RI30)
      SXY = Y * DX*(Y*RI30 - Y0 *RI31)
      SYX = Y * DX*(Y*RI31 - Y0 *RI30)
      SYY = Y *    (  RI11 + YY2*RI31-YYP*(RI30+RI32))
 
c------
c build the stress tensor
c------

      if(Iopt.ne.1) then

        QXXX = - Y6DX3 * RI50
        QXXY = - Y6DX2 * (Y*RI50 - Y0*RI51)
        QXYX =   QXXY
        QXYY = - Y6DX  * (Y02*RI52+Y2*RI50 - 2.0*YYP*RI51)
        QYXX = - Y6DX2 * (Y  *RI51-Y0*RI50)
        QYXY = - Y6DX  * ((Y02+Y2)*RI51 - YYP*(RI52+RI50))
        QYYX = QYXY
        QYYY = - Y66 * ( Y3     *  RI51
     +                 - Y03    *  RI52
     +                 - Y*YYP  * (RI50+2.0*RI52)
     +                 + Y0*YYP * (RI53+2.0*RI51) )
         PXX = QYXX
         PXY = QYXY
         PYX = - Y6DX * (Y2*RI52 - 2.0D0 * YYP*RI51 + Y02*RI50)
         PYY = - Y66  * (  Y0*YYP * (2.0D0 * RI52+RI50)
     +                   - Y03    * RI51
     +                   - Y *YYP * (2.0*RI51+RI53)
     +                   + Y3     * RI52  )

      End If
c------

c-----
c done
c-----

      return
      end
