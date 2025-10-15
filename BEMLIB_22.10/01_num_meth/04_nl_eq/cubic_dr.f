      program cubic_dr

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
c Roots of the cubic equation:
c
c  x**3 + a*x**2 + b*x + c = 0
c
c where a,b,c are real coefficients.
c
c The roots are computed
c using Cardano's analytical formulas
c------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      Null = 0
      None = 1
      Ntwo = 2
      Nthr = 3

c------
c input
c------

  98  Continue

      write (6,*)
      write (6,*) " Roots of: x**3 + a*x**2 + b*x + c = 0"
      write (6,*)
      write (6,*) " Please enter the coefficients a, b, c"
      write (6,*) "        one by one; followed by Return"
      write (6,*) "        after each entry;"
      write (6,*) "        enter 99 for any one to quit"
      write (6,*) " -------------------------------------"

      read  (5,*) a
      If(a.eq.99) Go to 99
      read  (5,*) b
      If(b.eq.99) Go to 99
      read  (5,*) c
      If(c.eq.99) Go to 99

c-------
c launch
c-------

      call cubic 
     +
     +  (a,b,c
     +  ,D
     +  ,x1,x2,x3
     +  ,prtr,prti
     +  )

c--------------------------------------
c one real, two complex conjugate roots
c--------------------------------------

      If(D.ge.0) then 

       write (6,*)
       write (6,*) "One real, two complex conjugate roots"
       write (6,*)
       write (6,*) 'Real and imaginary part of the roots:'
       write (6,*)

       write (6,101) None,x1,Null
       write (6,101) Ntwo,prtr, prti
       write (6,101) Nthr,prtr,-prti

c---
c error
c---

       error1r = x1**3+a*x1**2+b*x1+c
       error1i = 0.0

       error2r = prtr**3 - 3.0*prtr*prti**2
     +          + a*(prtr**2-prti**2)
     +          + b*prtr
     +          + c
       error2i =   3.0*prtr**2 * prti - prti**3
     +           + 2.0*a*prtr*prti
     +           + b*prti

       error3r =  error2r
       error3i = -error2i

       write (6,*)
       write (6,*) 'Real and imaginary part of the error:'
       write (6,*)

       write (6,101) None,error1r,error1i
       write (6,101) Ntwo,error2r,error2i
       write (6,101) Nthr,error3r,error3i

c------------------
c  three real roots
c------------------

      Else       ! three real roots

       write (6,*)
       write (6,*) "Three real roots:"
       write (6,*)

       write (6,101) None,x1
       write (6,101) Ntwo,x2
       write (6,101) Nthr,x3

c---
c error
c---

       error1 = x1**3+a*x1**2+b*x1+c
       error2 = x2**3+a*x2**2+b*x2+c
       error3 = x3**3+a*x3**2+b*x3+c

       write (6,*)
       write (6,*) "Error:"
       write (6,*)

       write (6,101) None,error1
       write (6,101) Ntwo,error2
       write (6,101) Nthr,error3

c-----------
      End If
c-----------

      Go to 98   ! Repeat

c-----
c Done
c-----

 99   Continue

 101  Format (1x,i1,2(1x,f15.10))

      Stop
      End
