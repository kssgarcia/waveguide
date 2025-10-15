      program quadc_dr

c=========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c========================================

c------------------------------------
c Roots of the quadratic equation:
c
c  a*x**2 + b*x + c = 0
c
c where a,b,c are complex coefficients
c
c The roots are computed
c using the quadratic formula
c------------------------------------

      Implicit Double Precision (a-h,o-z)

c------
c constants
c------

      None = 1
      Ntwo = 2

c------
c input
c------

  98  Continue

      write (6,*)
      write (6,*) " Roots of: a*x**2 + b*x + c = 0"
      write (6,*)
      write (6,*) " Please enter the coefficients a, b, c"
      write (6,*) "        one by one; followed by Return"
      write (6,*) "        after each entry;"
      write (6,*) "        enter 99 for any one to quit"
      write (6,*) " -------------------------------------"

      read  (5,*) ar,ai
      if(ar.eq.99) Go to 99
      if(ai.eq.99) Go to 99
      read  (5,*) br,bi
      if(br.eq.99) Go to 99
      if(bi.eq.99) Go to 99
      read  (5,*) cr,ci
      if(cr.eq.99) Go to 99
      if(ci.eq.99) Go to 99

c-------
c launch
c-------

      call quadc
     +
     +  (ar,ai,br,bi,cr,ci
     +  ,ZR1,ZI1
     +  ,ZR2,ZI2
     +  )


c--------------------------------------
c one real, two complex conjugate roots
c--------------------------------------

       write (6,101) None,ZR1,ZI1
       write (6,101) Ntwo,ZR2,ZI2

c-----
c Done
c-----

 99   Continue

 101  Format (1x,i1,2(1x,f15.10))

      Stop
      End
