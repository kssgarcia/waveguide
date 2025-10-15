      subroutine jacobi (a,b,n,t,rjac)

c==================================
c FDLIB, BEMLIB
c
c evaluation of a Jacobi polynomial
c
c jac = J^(a,b)_n(t)
c
c by recursion
c==================================

      Implicit Double Precision (a-h,o-z)

      Dimension p(100)

c---
c use the recursion formula
c---

      rjac=1.0;

      if(n.gt.0) then

        al = 0.5D0*(a+b+2)
        be= (b-a)/(a+b+2)
        p(1)= al*(t-be)

         if(n.gt.1) then

          al = (2+a+b+1)*(2+a+b+2)/(4*(a+b+2) )
          be = (b**2-a**2)/((2+a+b)*(2+a+b+2))
          ga = (1+a)*(1+b)*(2+a+b+2)/(2*(a+b+2)*(2+a+b))

         p(2)= al*(t-be)*p(1) - ga;

         if(n.gt.2) then

         Do i=2,n-1
           al = (2*i+a+b+1)*(2*i+a+b+2)/(2*(i+1)*(i+a+b+1) )
           be = (b**2-a**2)/((2*i+a+b)*(2*i+a+b+2))
           ga = (i+a)*(i+b)*(2*i+a+b+2)/((i+1)*(i+a+b+1)*(2*i+a+b))
           p(i+1)= al*(t-be)*p(i) - ga*p(i-1)
         End Do 

         end if
         end if

         rjac=p(n)

      end if

c=== Lobatto
c Lo = 0.5*(n+2)*jac
c===
c Lo1 = 3*t
c Lo5 = 1/8  * (693*t^4-630*t^2+105)*t
c Lo6 = 1/16 * (3003*t^6-3465*t^4+945*t^2-35)
c
c=== Legendre
c L = rjac
c L6 =  1/16 * (231*t^6 - 315*t^4 + 105*t^2 - 5)
c=====

c-----
c done
c-----

      return
      end
