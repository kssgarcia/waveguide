      subroutine deflation 
     +
     +  (Idfl
     +  ,Isym
     +  ,U,V
     +  ,CF1,CF2,CF3,CF4
     +  ,xcnt,ycnt
     +  )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c----------------------------------------
c Compute deflation integrals
c
c SYMBOLS:
c -------
c
c ARL: arc length
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension X(0:900),Y(0:900)
      Dimension U(0:900),V(0:900)

      Dimension  XC(900), YC(900),  R(900),   S(900)
      Dimension TH1(900),TH2(900),TH3(900),ORNT(900)

      Dimension ZZ(20),WW(20)

c---
c common blocks
c---

      common/XXYY/X,Y
      common/ARCC/XC,YC,R,S,TH1,TH2,TH3,ORNT

      common/ancI1/NSG,NSG1,NSG2,NSGM,NSGQ
      common/ancI3/NGL

      common/ZZWW/ZZ,WW
      common/piii/pi,pih,piq,pi2,pi4,pi6,pi8,srpi

c------------------------------
c reflect a symmetric interface
c------------------------------

      if(Isym.eq.2) then

        Do i=1,NSGQ
         U(NSGM-i+1) = -U(i)
         V(NSGM-i+1) =  V(i)
         U(NSGM+i-1) = -U(i)
         V(NSGM+i-1) = -V(i)
         U(NSG -i+2) =  U(i)
         V(NSG -i+2) = -V(i)
        End Do

      end if

c-----------
c initialize
c-----------

      ARL  = 0.0D0          ! arc length
      PRJ1 = 0.0D0          ! v-dot-n

      if(Idfl.eq.2) then    ! complete deflation
        PRJ2 = 0.0D0        ! int(u)
        PRJ3 = 0.0D0        ! int(v)
        PRJ4 = 0.0D0
        PRJ5 = 0.0D0
        PRJ6 = 0.0D0
        PRJ7 = 0.0D0
        PRJ8 = 0.0D0
      end if

c--------------------
c deflation integrals
c--------------------

      Do J=1,NSG

        J1 = J+1

        ST1 = 0.5D0*(TH3(J) +TH2(J))
        DT1 = 0.5D0*(TH3(J) -TH2(J))
        ST2 = 0.5D0*(TH2(J1)+TH1(J1))
        DT2 = 0.5D0*(TH2(J1)-TH1(J1))

        SUD = 0.5D0*(U(J1)+ U(J))
        DUD = 0.5D0*(U(J1)- U(J))
        SVD = 0.5D0*(V(J1)+ V(J))
        DVD = 0.5D0*(V(J1)- V(J))

        SUM1 = 0.0D0

        if(Idfl.eq.2) then
          SUM2 = 0.0D0
          SUM3 = 0.0D0
          SUM4 = 0.0D0
          SUM5 = 0.0D0
          SUM6 = 0.0D0
          SUM7 = 0.0D0
          SUM8 = 0.0D0
        end if

c---
c Gauss--Legendre  quadrature
c---

        Do k=1,NGL

          UU = SUD + DUD*ZZ(k)
          VV = SVD + DVD*ZZ(k)

          AN1 = ST1+DT1*ZZ(k)
          CS1 = Dcos(AN1)
          SN1 = Dsin(AN1)
          XX1 = XC(J)+R(J)*CS1
          YY1 = YC(J)+R(J)*SN1
          vx1 = CS1*ORNT(J)
          vy1 = SN1*ORNT(J)

          AN2 = ST2+DT2*ZZ(k)
          CS2 = Dcos(AN2)
          SN2 = Dsin(AN2)
          XX2 = XC(J1) + R(J1)*CS2
          YY2 = YC(J1) + R(J1)*SN2
          vx2 = CS2*ORNT(J1)
          vy2 = SN2*ORNT(J1)

          XX = 0.5D0*(XX1+XX2)
          YY = 0.5D0*(YY1+YY2)
          vx = 0.5D0*(vx1+vx2)
          vy = 0.5D0*(vy1+vy2)

          www = WW(k)

          SUM1 = SUM1 +(UU*vx+VV*vy)*www

          if(Idfl.eq.2) then
            SUM2 = SUM2 + UU           * www
            SUM3 = SUM3 + VV           * www
            SUM4 = SUM4 + XX**2        * www
            SUM5 = SUM5 + YY**2        * www
            SUM6 = SUM6 + XX           * www
            SUM7 = SUM7 + YY           * www
            SUM8 = SUM8 +(-YY*UU+XX*VV)* www
          end if

        End Do

        cf = 0.5D0*(DT1*R(J)*ORNT(J)+DT2*R(J1)*ORNT(J1))

        ARL  = ARL  +      cf
        PRJ1 = PRJ1 + SUM1*cf

        if(Idfl.eq.2) then
          PRJ2 = PRJ2 + SUM2*cf
          PRJ3 = PRJ3 + SUM3*cf
          PRJ4 = PRJ4 + SUM4*cf
          PRJ5 = PRJ5 + SUM5*cf
          PRJ6 = PRJ6 + SUM6*cf
          PRJ7 = PRJ7 + SUM7*cf
          PRJ8 = PRJ8 + SUM8*cf
        end if

      End Do

c----------------
c post processing
c----------------

      ARL = 2.0D0*ARL         ! factor of 2 bec mult by DTH in cf
      CF1 = PRJ1/ARL

      if(Idfl.eq.2) then

        CF2  = PRJ2/ARL
        CF3  = PRJ3/ARL
        xcnt = PRJ6/ARL
        ycnt = PRJ7/ARL
        PRJ8 = PRJ8 + ycnt*PRJ2 - xcnt*PRJ3
        Rot_norm = PRJ4+PRJ5
     +           -2.0*(xcnt*PRJ6+ycnt*PRJ7)
     +           +(xcnt**2+ycnt**2)*ARL
        CF4  = PRJ8/Rot_norm

      end if

c-----
c done
c-----

      return
      end
