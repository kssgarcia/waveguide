      function leg_poly(n,x,p)

c===========================================
c FDLIB
c
c Evaluate the Legendre polynomial Pn(x)
c by recursion
c===========================================

      Implicit Double Precision (a-h,o-z)

      if(n.eq.0) then
       p = 1.0D0
      else if(n.eq.1) then
       p = x
      else
       p0 = 1.0D0
       p1 = x
       Do i=1,n-1
        p2 = (2.0D0*i+1.0D0)/(i+1.0D0) *x*p1 - i/(i+1.0D0) * p0
        p0 = p1
        p1 = p2
       End Do
       p = p2
      end if

c---
c done
c---

      return
      end
