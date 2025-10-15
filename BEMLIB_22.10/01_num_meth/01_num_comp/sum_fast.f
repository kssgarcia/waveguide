      program fast_sum

c==========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c===============================================
c This program accompanies the book:
c           C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c        Oxford University Press
c===============================================

c----------------------------------------------
c Fast summation of an infinite sum: sum{s(i)}
c whose terms decay like
c
c  s(i) ~ 1/i**2
c
c SYMBOLS
c -------
c
c kmax: number of terms in the computed x sequence
c       (set to 7)
c-------------------------------------------------

c------------------------------------------------
c          C. Pozrikidis
c Numerical Computation in Science and Engineering
c      Oxford University Press
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension M(100),x(100),a(100),t(100)

      Parameter(kmax=7)

c---
c constants
c---

      pi = 3.14159265358 D0

      None = 1

c---
c menu
c---

  97  Continue

      write (6,*)
      write (6,*) " Choose the sequence:"
      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 for s(i)=1/i**2"
      write (6,*) " 2 for s(i)=1/(i*(i+2))"
      write (6,*) " 3 for s(i)=i**2/(1+i**4)"
      write (6,*) " 4 for s(i)=1/(2*i-1)**2"
      write (6,*) " 5 for s(i)=1/i**s"
      write (6,*) " 0 to quit"
      write (6,*) "----------"
      read  (5,*) menu

      if(menu.eq.0) Go to 99

c--------------------------------
c print the exact value, if known
c--------------------------------

      if(menu.eq.2) then
         write (6,*) 
         write (6,*) " The exact value is 0.75"
         write (6,*) 
      elseif(menu.eq.4) then
         exact = pi**2/8.0D0
         write (6,*) 
         write (6,*) " The exact value is ",exact
         write (6,*) 
      elseif(menu.eq.5) then
         write (6,*) " Enter the exponent s"
         read  (5,*) s
      end if

c----------
c inquiries
c----------

      write (6,*)
      write (6,*) " Enter the first integer factor N (try 8)"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) N

      if(N.eq.0) Go to 99

      write (6,*)
      write (6,*) " Enter the integer p (try 2)"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) Np

      if(Np.eq.0) Go to 99

c----------------------------
c max term of the x(k) series
c----------------------------

      M(1) = N

      Do k=2,kmax
       M(k) = M(k-1)*Np
      End Do

c------------------------
c compute the x(k) series
c------------------------

      sum = 0.0D0

      Do i=1,M(kmax)

        if(menu.eq.1) then
          term = 1.0D0/i**2
        elseif(menu.eq.2) then
          term = 1.0D0/(i*(i+2.0D0))
        elseif(menu.eq.3) then
          is = i**2
          term = 1.0D0/(1.0D0/is+is)
        elseif(menu.eq.4) then
          term = 1.0D0/(2.0D0*i-1.0D0)**2
        elseif(menu.eq.5) then
          term = 1.0D0/i**s
        end if

        sum = sum+term

        if(i.eq.M(1)) then
           t(1) = term
           x(1) = sum - 0.5*term
        elseif(i.eq.M(2)) then
           t(2) = term
           x(2) = sum - 0.5*term
        elseif(i.eq.M(3)) then
           t(3) = term
           x(3) = sum - 0.5*term
        elseif(i.eq.M(4)) then
           t(4) = term
           x(4) = sum - 0.5*term
        elseif(i.eq.M(5)) then
           t(5) = term
           x(5) = sum - 0.5*term
        elseif(i.eq.M(6)) then
           t(6) = term
           x(6) = sum - 0.5*term
        elseif(i.eq.M(7)) then
           t(7) = term
           x(7) = sum - 0.5*term
        endif
      End Do

c---------------------
c Aitken extrapolation
c---------------------

      Do k=2,kmax-1
        a(k) = (x(k-1)*x(k+1)-x(k)**2)
     +       / (x(k+1)-2.0D0*x(k)+x(k-1))
c    +       + 0.5*t(k)
      End Do

c----------------------------------------
c print the raw and extrapolated sequence
c----------------------------------------

      write (6,101) None,M(1),x(1)

      Do k=2,kmax-1
       x(k) = x(k)+0.5D0*t(k)
       write (6,101) k,M(k),x(k),a(k)
      End Do

      write (6,101) kmax,M(kmax),x(kmax)

      Go to 97

c-----
c Done
c-----

 99   Continue

 100  Format(4(1x,f15.10))
 101  Format(1x,i3,1x,i9,4(1x,f15.10))

      stop
      end
