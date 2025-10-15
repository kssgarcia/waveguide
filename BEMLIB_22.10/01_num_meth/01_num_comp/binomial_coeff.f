      program binomial_coeff

c====================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c====================================

c------------------------------------------------
c This program accompanies the book:
c             C. Pozrikidis
c "Numerical Computation in Science and Engineering"
c        Oxford University Press
c------------------------------------------------

c----------------------------------------------------
c Computation of the binomial coefficient in two ways
c----------------------------------------------------

 98   Continue

      write (6,*)
      write (6,*) " Will compute the (n/m) combinatorial"
      write (6,*)

 97   Continue

      write (6,*) " Please enter n and m (n>=m)"
      write (6,*) " -1 for either one to quit"
      write (6,*) "-------------------------"
      read (5,*) n,m

      if(n.eq.-1) Go to 99
      if(m.eq.-1) Go to 99

c---
c trap
c---

      if(m.gt.n) then
       write (6,*)
       write (6,*) "Invalid input; please try again"
       write (6,*) "m cannot be greater than n"
       write (6,*)
       Go to 97
      end if

c-------------
c direct method
c--------------

      nf = 1      ! compute n!

      Do i=1,n
       nf = nf*i
      End Do

      mf = 1      ! compute m!

      Do i=1,m
       mf = mf*i
      End Do

      nmmf = 1    ! compute (n-m)!

      Do i=1,n-m
       nmmf = nmmf*i
      End Do

      nmc = nf/(mf*nmmf)

      write (6,*) 
      write (6,*) " Binomail coefficient by the direct method:"
      write (6,*) nmc

c----------------------------------------------
c  better method according to algorithm (1.4.7)
c----------------------------------------------

c---
c find the minimum of m and n-m
c---

      l = m
      if((n-m).lt.l) l = n-m

c---
c some magic
c---

      nmc = 1
      Do k=1,l
       nmc = nmc*(n-k+1)/k
      End Do

      write (6,*) 
      write (6,*) " Binomial coefficient by the improved method:"
      write (6,*) nmc

      Go to 98

c-----
c done
c-----

  99  Continue

      stop
      end
