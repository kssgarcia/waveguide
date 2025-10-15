      subroutine crv_max (crvmax)

c======================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c======================================

c----------------------------------------
c Computate the maximum curvature
c
c Isym = 1 indicates that the line is symmetric
c          with respect to the mid-point
c          In that case, minimum curvature
c          is assumed to occur at the mid-point
c
c          Modify as necessary
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension X(0:900),Y(0:900)

      Dimension  XC(900), YC(900),  R(900),   S(900)
      Dimension TH1(900),TH2(900),TH3(900),ORNT(900)

      common/ancI1/NSG,NSG1,NSG2,NSGM,NSGQ
      common/XXYY/X,Y
      common/ARCC/XC,YC,R,S,TH1,TH2,TH3,ORNT

c-----------------------------------------
c search for minimum curvature
c followed by quadratic interpolation
c-----------------------------------------

      F = 100000.0
      M = 0

      Do i=1,NSG
        tmp  = abs(R(i))
        if(tmp.lt.F) then
           F = tmp
           M = i
        end if
      End Do

      Ma = M-1
      M1 = M
      M2 = M+1

c----------------------------------
c compute interpolation arc lengths
c----------------------------------

      xx0 = - Dsqrt((X(Ma)-X(M1))**2+(Y(Ma)-Y(M1))**2)
      xx1 = 0.0
      xx2 =   Dsqrt((X(M2)-X(M1))**2+(Y(M2)-Y(M1))**2)

      yy0 = Dabs(R(Ma))
      yy1 = F
      yy2 = Dabs(R(M2))

      tm2 = (yy2-yy1)/(xx2-xx1)
      tm1 = (yy0-yy1)/(xx0-xx1)
      bbb = (tm2-tm1)/(xx2-xx0)

      ccc = tm2-bbb*(xx2-xx1)

      Rmin = - 0.25*ccc**2/bbb + yy1   ! minimum radius of curvature

      crvmax = 1.0/Rmin

c-----
c done
c-----

      return
      end
