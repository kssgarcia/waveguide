      subroutine cubic
     +
     +  (a,b,c
     +  ,D
     +  ,x1,x2,x3
     +  ,prtr,prti
     +  )

c=================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=================================

c------------------------------------
c Roots of the cubic equation:
c
c  x^3 + a*x^2 + b*x + c = 0
c
c where a, b, c are real coefficients.
c
c The roots are computed
c using Cardano's analytical formulae
c
c See:
c
c           John W Harris and Horst Stocker,
c ``Handbook of Mathematics and Computational Science''
c                Springer (1998).
c
c john.harris@yale.edu
c------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi = 3.14159265358D0

      oot   = 1.0D0/3.0D0
      opf   = 1.5D0
      three = 3.0D0
      srth  = Dsqrt(three)

c--------
c prepare
c--------

      p = b-a*a/3.0D0
      q = c+2.0D0*a*a*a/27.0D0-a*b/3.0D0
      D = (p/3.0D0)**3+(q/2.0D0)**2        ! discriminant

c-------------------------------------------
c one real root, two complex conjugate roots
c-------------------------------------------

      if(D.ge.0) then 

       srd  = Dsqrt(D)
       tmp  =  -0.5D0*q+srd
       u    =  Dabs(tmp)**oot
       if(tmp.lt.0) u = -u
        tmp  =  -0.5D0*q-srd
        v    =  Dabs(tmp)**oot
       if(tmp.lt.0) v = -v

       x1 = -a/3.0D0+u+v

       prtr = -a/3.0D0-0.5D0*(u+v)
       prti = srth*0.5D0*(u-v)

c-----------------
c three real roots
c-----------------

      else

       cosphi = -0.5D0*q/(Dabs(p)/3.0D0)**opf
       phi = Dacos(cosphi)

       cf = 2.0D0*sqrt(Dabs(p)/3.0D0)

       x1 = -a/3.0D0 + cf*Dcos(phi/3.0D0)
       x2 = -a/3.0D0 - cf*Dcos((phi-pi)/3.0D0)
       x3 = -a/3.0D0 - cf*Dcos((phi+pi)/3.0D0)

c-----------
      end if
c-----------

c-----
c done
c-----

      return
      end
