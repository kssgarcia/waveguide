      subroutine sgf_3x
     +
     +  (X,Y
     +  ,X0,Y0
     +  ,QRxx1,QRxs1,QRsx1
     +  ,QRss1,QRff1
     +  ,QIxf1
     +  ,QIfx1,QIsf1,QIfs1
     +  )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c--------------------------------------
c Green's function of three-dimensional 
c Stokes flow in an axisymmetric domain
c for use with the boundary-integral
c equation
c-------------------------------------

      Implicit Double Precision (a-h,o-z)

c--------
c prepare
c--------

      Y2  =  Y*Y
      Y3  =  Y*Y2
      Y4  = Y2*Y2

      Y02 =  Y0*Y0
      Y03 =  Y0*Y02
      Y04 = Y02*Y02

      SY  = Y+Y0
      SY2 = SY*SY
      YYS = Y2+Y02

      YY2 = Y*Y0
      YYR = Dsqrt(YY2)
      YY3 = YY2*YYR
      YY5 = YYR*YYR*YY3

      DX  = X-X0
      DX2 = DX*DX
      DX4 = DX2*DX2

      DY  = Y-Y0
      DY2 = DY*DY

      DR2 = DX2+DY2
      DR  = Dsqrt(DR2)

c      DXYY = DX2+YYS

c------------------------
c additional preparations
c------------------------

      RK2  = 4.0D0*YY2/(DX2+SY2)
      RK   = Dsqrt(RK2)
      RK3  = RK2*RK
      RK4  = RK2*RK2
      RK6  = RK4*RK2

      RK2P = 1.0D0-RK2
      RK2H = 1.0D0-0.50D0*RK2

      call ell_int (RK,F,E)

      RL10 = F
      RL11 = 2.0D0*( RK2H*F-E ) / RK2
      RL12 = ( (3.0*RK4-8.0*RK2+8.0)*F-8.0*RK2H*E)/(3.0*RK4)

      FCTR = 2.0*RK/YYR

      RI10 = FCTR*RL10
      RI11 = FCTR*RL11
      RI12 = FCTR*RL12

c      RI11 = RK*(DXYY*F-(DX2+SY2)*E)/YY3

      RL30 = E/RK2P
      RL31 = 2.0*(-F+ RK2H/RK2P * E ) / RK2
      RL32 = 4.0*(RK2H**2 * RL30 - 2.0*RK2H*RL10 + E ) / RK4
      RL33 = 8.0*(RK2H**3 * RL30 - 3.0*RK2H**2 * F
     +      + 3.0*RK2H    * E - (4.0*RK2H*E-RK2P*F)/3.0 ) / RK6

      FCTR = RK3/(2.0*YY3)

      RI30 = FCTR*RL30
      RI31 = FCTR*RL31
      RI32 = FCTR*RL32
      RI33 = FCTR*RL33

c------------
c      RI31 = RK*(-F+DXYY*E/DR2)/YY3
c      RI32 = RK*(-DXYY*F+(DX4+2.0*DX2*YYS+Y4+Y04)*E/DR2)/YY5
c------------

      QRxx1 = Y*   (  RI11 + DX2 * RI31)
      QRxs1 = Y*DX*(Y*RI31 - Y0  * RI32)
      QRsx1 = Y*DX*(Y*RI32 - Y0  * RI31)
      QRss1 = Y*(RI12 + YYS*RI32 - YY2*(RI33+RI31))
      QRff1 = Y*(RI12            + YY2*(RI31-RI33))

      QIxf1 = - Y*DX*Y0*(RI32-RI30)
      QIfx1 =   Y*DX*Y *(RI30-RI32)
      QIsf1 =   Y*(-RI10+RI12-Y02*(RI30-RI32)+YY2*(RI31-RI33))
      QIfs1 =   Y*( RI10-RI12+Y2 *(RI30-RI32)-YY2*(RI31-RI33))
    
c-----
c Done
c-----

 100  Format (3(1x,f15.10))

      Return
      End
