      subroutine caps_2d_rhe
     +
     +    (strxx
     +    ,strxy
     +    ,stryy
     +    )

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c-----------------------------------------------
c Computes the effective stress tensor 
c of an infinite dilute suspension
c
c (FX,FY) is generated in subroutine caps_2d_slp
c-----------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  X(0:900), Y(0:900)
      Dimension  U(0:900), V(0:900)
      Dimension FX(0:900),FY(0:900)

      Dimension  XC(900), YC(900),  R(900),   S(900)
      Dimension TH1(900),TH2(900),TH3(900),ORNT(900)

      Dimension ZZ(20),WW(20)

c--------------
c common blocks
c--------------

      common/XXYY/X,Y
      common/ARCC/XC,YC,R,S,TH1,TH2,TH3,ORNT
      common/UUVV/U,V

      common/DDFF/FX,FY

      common/ancR2/vs1,vs2,vsr,vsr1,vsr2,rkap,rkap2,eps
      common/ancI1/NSG,NSG1,NSG2,NSGM,NSGQ
      common/ancI3/NGL
      common/ZZWW/ZZ,WW

c--------
c prepare
c--------

      vsdiff = vs1-vs2

c-----------
c initialize
c-----------

      strxx = 0.0D0
      strxy = 0.0D0
      stryx = 0.0D0
      stryy = 0.0D0

      arl = 0.0D0

c---
c integrate over the blended arc
c---


      Do i=1,NSG

      i1 = i+1

      sxx = 0.0
      sxy = 0.0
      syx = 0.0
      syy = 0.0

      XC1 =   XC(i)
      YC1 =   YC(i)
      RD1 =    R(i)
      T11 =  TH2(i)
      T21 =  TH3(i)
      OR1 = ORNT(i)
      FX1 =   FX(i)
      FY1 =   FY(i)
      UX1 =    U(i)
      UY1 =    V(i)

      XC2 =   XC(i1)
      YC2 =   YC(i1)
      RD2 =    R(i1)
      T12 =  TH1(i1)
      T22 =  TH2(i1)
      OR2 = ORNT(i1)
      FX2 =   FX(i1)
      FY2 =   FY(i1)
      UX2 =    U(i1)
      UY2 =    V(i1)

      ST1 = 0.5D0*(T21+T11)
      DT1 = 0.5D0*(T21-T11)

      ST2 = 0.5D0*(T22+T12)
      DT2 = 0.5D0*(T22-T12)

      SFX = 0.5D0*(FX2+FX1)
      DFX = 0.5D0*(FX2-FX1)

      SFY = 0.5D0*(FY2+FY1)
      DFY = 0.5D0*(FY2-FY1)

      SUX = 0.5D0*(UX2+UX1)
      DUX = 0.5D0*(UX2-UX1)

      SUY = 0.5D0*(UY2+UY1)
      DUY = 0.5D0*(UY2-UY1)

      Do k=1,NGL

      AN1 = ST1 + DT1*ZZ(k)
      CS1 = Dcos(AN1)
      SN1 = Dsin(AN1)
      XX1 = XC1 + RD1 * CS1
      YY1 = YC1 + RD1 * SN1

      AN2 = ST2 + DT2*ZZ(k)
      CS2 = Dcos(AN2)
      SN2 = Dsin(AN2)
      XX2 = XC2 + RD2 * CS2
      YY2 = YC2 + RD2 * SN2

      XX = 0.5*(XX1+XX2)
      YY = 0.5*(YY1+YY2)

      FXX = SFX + DFX*ZZ(k)
      FYY = SFY + DFY*ZZ(k)

      UX = SUX + DUX*ZZ(k)
      UY = SUY + DUY*ZZ(k)

      vnx = 0.5D0*(CS1*OR1 + CS2*OR2)
      vny = 0.5D0*(SN1*OR1 + SN2*OR2)
      par = sqrt(vnx**2+vny**2)
      vnx = vnx/par
      vny = vny/par

      sxx = sxx + (fxx*xx - vsdiff*(ux*vnx+ux*vnx))*WW(k)
      sxy = sxy + (fxx*yy - vsdiff*(ux*vny+uy*vnx))*WW(k)
      syx = syx + (fyy*xx - vsdiff*(uy*vnx+ux*vny))*WW(k)
      syy = syy + (fyy*yy - vsdiff*(uy*vny+uy*vny))*WW(k)

      End Do

      cf = 0.5D0*( DT1*OR1*RD1 + DT2*OR2*RD2 )

      strxx = strxx + sxx*cf
      strxy = strxy + sxy*cf
      stryx = stryx + syx*cf
      stryy = stryy + syy*cf

      arl = arl+cf

      End Do

      arl = 2.0D0*arl

c---------
c printing
c---------

c     write (6,*)
c     write (6,*) " Effective stresses"
c     write (6,*)
c     write (6,100) strxx,strxy
c     write (6,100) stryx,stryy
c     write (6,*)
c     write (6,*) " Arc length"
c     write (6,*)
c     write (6,100) arl

c-----
c done
c-----

 100  Format (3(1x,f15.8))
 101  Format (1x,i3,5(1x,f15.8))

      return
      end

c=====================================

      subroutine channel 
     +                   (IQPD
     +                   ,RL
     +                   ,prgr
     +                   ,flrt
     +                   )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement.
c----------------------------------------

c----------------------------------------
c Compute the disturbance pressure drop 
c for channel flow
c
c (FX,FY) is generated in subroutine drop_2d_slp
c
c SYMBOLS:
c -------
c
c prgr:  pressure gradient
c
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  XC(900), YC(900),  R(900),   S(900)
      Dimension TH1(900),TH2(900),TH3(900),ORNT(900)

      Dimension FX(0:900),FY(0:900)

      Dimension ZZ(20),WW(20)

c---
c common blocks
c---

      common/ARCC/XC,YC,R,S,TH1,TH2,TH3,ORNT
      common/DDFF/FX,FY
      common/hhhh/h,hh,h2,h3,h4,hs
      common/ancI1/NSG,NSG1,NSG2,NSGM,NSGQ
      common/ancI3/NGL
      common/ZZWW/ZZ,WW

c---
c evaluate the interfacial integral
c by quadrature
c---

      Bin = 0.0

      Do J=1,NSG

        J1 = J+1
        SFX = 0.5*(FX(J1)+FX(J))
        DFX = 0.5*(FX(J1)-FX(J))
        STH = 0.5*(TH3(J)+TH2(J))
        DTH = 0.5*(TH3(J)-TH2(J))

        Sum = 0.
        Do k = 1,NGL
          FXX = SFX + DFX*ZZ(k)
          TH  = STH + DTH*ZZ(k)
          SN  = SIN(TH)
          YY  = YC(J) + R(J)*SN
          Sum = Sum +FXX*(1.0-YY**2/hs)*WW(k)
        End Do

        cf = ORNT(J)*DTH*R(J)
        Bin = Bin + Sum*cf

      End Do

      If(IQPD.eq.0) then
        prgr = 0.
        flrt = -0.5*hs/RL * Bin
      Else
        prgr = 3.0/(4.0*RL*h) * Bin
        flrt = 0.
      End If

c-----
c Done
c-----

      Return
      End
