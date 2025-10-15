      subroutine em_2d_pdfl
     +
     +    (NSG
     +    ,IQPD
     +    ,RL
     +    ,prgr
     +    ,flrt
     +    )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c----------------------------------------------
c Compute the disturbance pressure drop
c of the distrurbance flow rate
c for channel flow
c
c (FX,FY) is generated in subroutine em_2d_slp
c
c SYMBOLS:
c -------
c
c prgr:  pressure gradient
c---------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension XC (200),YC (200),R  (200),S   (200)
      Dimension TH1(200),TH2(200),TH3(200),ORNT(200)

      Dimension FX(0:200),FY(0:200)

      Dimension ZZ(20),WW(20)

c---
c common blocks
c---

      common/ARCC/XC,YC,R,S,TH1,TH2,TH3,ORNT
      common/DDFF/FX,FY
      common/hhhh/h,hh,h2,h3,h4,hs
      common/ancI3/NGL
      common/ZZWW/ZZ,WW

c---
c evaluate an interfacial integral
c by quadrature
c---

      Bin = 0.0D0

      Do J=1,NSG

        J1 = J+1
        SFX  = 0.5D0*(FX(J1)+FX(J))
        DFX  = 0.5D0*(FX(J1)-FX(J))
        STH  = 0.5D0*(TH3(J)+TH2(J))
        DTH  = 0.5D0*(TH3(J)-TH2(J))

        Sum  = 0.0D0

        Do k=1,NGL
          FXX = SFX + DFX*ZZ(k)
          TH  = STH + DTH*ZZ(k)
          SN  = Dsin(TH)
          YY  = YC(J) + R(J)*SN
          Sum = Sum +FXX*(1.0-YY**2/hs)*WW(k)
        End Do

        cf = ORNT(J)*DTH*R(J)
        Bin = Bin + Sum*cf

      End Do

      if(IQPD.eq.0) then
        prgr = 0.0D0
        flrt = -0.5D0*hs/RL * Bin
      else
        prgr = 3.0D0/(4.0D0*RL*h) * Bin
        flrt   = 0.0D0
      end if

c-----
c done
c-----

      return
      end
