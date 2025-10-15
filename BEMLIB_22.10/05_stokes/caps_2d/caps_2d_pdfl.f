      subroutine caps_2d_pdfl
     +
     +    (IQPD
     +    ,RL
     +    ,prgr
     +    ,flrt
     +    )

c===========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licencing agreement
c===========================================

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
c evaluate an interfacial integral
c by quadrature
c---

      Bin = 0.0D0

      Do J=1,NSG

        J1 = J+1
        SFX = 0.5D0*(FX(J1)+FX(J))
        DFX = 0.5D0*(FX(J1)-FX(J))
        STH = 0.5D0*(TH3(J)+TH2(J))
        DTH = 0.5D0*(TH3(J)-TH2(J))

        Sum = 0.0D0

        Do k=1,NGL
          FXX = SFX + DFX*ZZ(k)
          TH  = STH + DTH*ZZ(k)
          SN  = Dsin(TH)
          YY  = YC(J) + R(J)*SN
          Sum = Sum +FXX*(1.0D0-YY**2/hs)*WW(k)
        End Do

        cf = ORNT(J)*DTH*R(J)
        Bin = Bin + Sum*cf

      End Do

      if(IQPD.eq.0) then
        prgr = 0.0D0
        flrt = -0.5D0*hs/RL * Bin
      else
        prgr = 3.0D0/(4.0D0*RL*h) * Bin
        flrt = 0.0D0
      end if

c-----
c done
c-----

      return
      end
