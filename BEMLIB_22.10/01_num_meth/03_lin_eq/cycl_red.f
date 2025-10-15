      subroutine cycl_red (n,p,a,b,c,s,x)

c========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c------------------------------------------------
c This program accompanies the book:
c
c          C. Pozrikidis
c Numerical Computation in Science and Engineering
c       Oxford University Press, 1998
c
c Cyclic reduction of a tridiagonal system:
c
c  T x = s 
c
c with constant diagonals
c
c Algorithm 3.4.2
c
c Number of equtions: n=2**p
c p is an integer
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension s(200),x(200),sig(200),snew(200)
      Integer p

c---
c prepare
c---

      al    = a
      be    = b
      ga    = c
      alast = a

      Nmax = n         ! size of descendant system

c----------------
c save s into sig
c----------------

      Do i=1,Nmax
       sig(i) = s(i)
      End Do

c----------
c reduction
c----------

      Do m=1,p

        anew = al - 2.0*be*ga/al
        bnew =      -   be**2/al
        cnew =      -   ga**2/al
        alast = alast - be*ga/al

        Do i=2,Nmax-2,2
          snew(i/2)  = -ga/al *sig(i-1)+sig(i)-be/al*sig(i+1)
        End Do

        snew(Nmax/2) = -ga/al*sig(Nmax-1) + sig(Nmax)

        al = anew
        be = bnew
        ga = cnew

        Nmax = Nmax/2

        write (6,100)
        write (6,100) Nmax

        Do i = 1,Nmax
          sig(i) = snew(i)
          write (6,100) i,sig(i)
        End Do

        write (6,100)

      End Do

c---------------------------
c solve for the last unknown
c---------------------------

      x(n) = sig(1)/alast

c------------------------------
c Provisional back substitution
c------------------------------

      x(n-1) = (s(n)-a*x(n))/c

      Do i=n-1,2,-1
       x(i-1) = (s(i)-a*x(i)-b*x(i+1))/c
      End Do
          
c-----
c Done
c-----

 100  Format (1x,i3,1x,f10.5)

      Return
      End
