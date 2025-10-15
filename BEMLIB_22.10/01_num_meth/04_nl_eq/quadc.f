      subroutine quadc
     +
     +  (AR,AI,BR,BI,CR,CI
     +  ,ZR1,ZI1
     +  ,ZR2,ZI2
     +  )

c==========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c-------------------------------
c  roots of a quadratic equation
c  with complex coefficients
c
c  A Z*Z + B Z + C = 0
c
c  quadratic formula:
c
c  Z = (-B+sqrt(B*B-4*A*C))/(2*A)
c-------------------------------

      Implicit Double Precision (a-h,o-z)

      Parameter(eps=0.000001)

c---
c linear equation
c---

      if((abs(AR).lt.eps).and.(abs(AI).lt.eps)) then
       tmp = BR*BR+BI*BI  ! inverse of B:
       BiR = BR/tmp
       BiI =-BI/tmp
       ZR1 = -(CR*BiR-CI*BiI)
       ZI1 = -(CR*BiI+CI*BiR)
       zR2 = 0.0D0
       zI2 = 0.0D0
       return
      end if

c---
c quadratic equation
c---

      DiscR = BR*BR-BI*BI-4.0D0*(AR*CR-AI*CI)
      DiscI = 2.0D0*BR*BI-4.0D0*(AR*CI+AI*CR)
      DiscM = Dsqrt(DiscR*DiscR+DiscI*DiscI)
      DiscA = acos(DiscR/DiscM)
      if(DiscI.lt.0) then
       DiscA = -DiscA
      end if
      tmp  = sqrt(DiscM)
      TMPR = tmp*Dcos(0.5D0*DiscA)
      TMPI = tmp*Dsin(0.5D0*DiscA)

      tmp = AR*AR+AI*AI     ! inverse of A
      AiR = AR/tmp
      AiI =-AI/tmp

      ZR1 = 0.5*(AiR*(-BR+TMPR)-AiI*(-BI+TMPI))
      ZI1 = 0.5*(AiI*(-BR+TMPR)+AiR*(-BI+TMPI))
 
      ZR2 = 0.5*(AiR*(-BR-TMPR)-AiI*(-BI-TMPI))
      ZI2 = 0.5*(AiI*(-BR-TMPR)+AiR*(-BI-TMPI))

c---
c done
c---

      return
      end
