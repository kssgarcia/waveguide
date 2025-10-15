      program leg_poly_dr

c===========================================
c FDLIB
c
c Evaluate the Legendre polynomial Pn(x)
c by recursion
c===========================================

      Implicit Double Precision (a-h,o-z)

      n = 4
      n = 3
      x = 0.9

      call leg_poly(n,x,p)

      p3 = (5*x**2-3)*x/2;
      p4 = (35*x**4-30*x**2+3)/8;

      write (6,100) x,p,p3,p4

c---
c product to sum
c---

      n = 3
      non = 1
      ntw = 2
      x = 0.45D0

      call leg_poly(non,x,p1)
      call leg_poly(ntw,x,p2)
      call leg_poly(n,x,pn)
      call leg_poly(n+1,x,pn1)
      call leg_poly(n-1,x,pna)
      call leg_poly(n+2,x,pn2)
      call leg_poly(n-2,x,pnb)
      prd1 = p1*pn
      prd2 = p2*pn

      c0 = (n+1.0D0)/(2.0D0*n+1.0D0)
      c1 = n/(2.0D0*n+1.0D0)

      c0r = (n+2.0D0)/(2.0D0*n+3.0D0)
      c1r = (n+1)/(2.0D0*n+3.0D0)

      c0s = n/(2.0D0*n-1.0D0)
      c1s = (n-1.0D0)/(2.0D0*n-1.0D0)


      prd1A = c0*pn1 + c1*pna

      c2 = 1.5D0*c0*c0r
      c0 = (3.0D0*c0*c1r+3.0D0*c1*c0s-1.0D0)/2.0D0
      c0 = n*(n+1.0D0)/(2*n-1.0D0)/(2*n+3.0D0)
      cb = 1.5D0*c1*c1s

      prd2A = c2*pn2 + c0*pn + cb*pnb

      write (6,100) x,prd1,prd1A
      write (6,100) x,prd2,prd2A
      

  100 format (1x,f10.5,1x,f10.5,1x,f10.5)

      stop
      end
