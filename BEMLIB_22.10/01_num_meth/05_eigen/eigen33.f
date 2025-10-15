      program eigen33

c=======================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c==========================================
c Eigenvalues of a 3x3 matrix with
c real elements
c
c The eigenvalues are computed by solving a cubic
c equation using Cardano's formula
c==========================================

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      Null = 0
      None = 1
      Ntwo = 2
      Nthr = 3

c----------------
c Read the matrix
c----------------

  98  Continue

      write (6,*) " Will compute the eigenvalues of a 3x3 matrix"
      write (6,*)
      write (6,*) " Please enter the matrix row-wise"
      write (6,*) " --------------------------------"

      read (5,*) a11,a12,a13
      read (5,*) a21,a22,a23
      read (5,*) a31,a32,a33

c     open (8,file="eigen33.dat")
c
c       Do i=1,3
c        read (8,*) a11,a12,a13
c        read (8,*) a21,a22,a23
c        read (8,*) a31,a32,a33
c       End Do
c
c      close (8)

c----------------------------------
c The characteristic polynomial will
c be placed in the form
c
c  x**3 + a*x**2 + b*x + c = 0
c
c where a,b,c are real coefficients.
c
c The roots will be computed
c using Cardano's analytical formulae
c
c compute the coefficients
c----------------------------------

      trace = a11+a22+a33

      Det = a11*(a22*a33-a23*a32)
     +     -a21*(a12*a33-a13*a32)
     +     +a31*(a12*a23-a13*a22)

      a = -trace
      b = (a11*a22-a21*a12)+(a22*a33-a23*a32)+(a11*a33-a13*a31)
      c = -Det

      call cubic
     +
     +   (a,b,c
     +   ,D
     +   ,x1,x2,x3
     +   ,prtr,prti
     +   )

c--------------------------------------
c one real, two complex conjugate roots
c--------------------------------------

      if(D.ge.0) then

       write (6,*)
       write (6,*) "One real, two complex conjugate roots"
       write (6,*)
       write (6,*) 'Real and imaginary part of the roots'
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
       write (6,*) 'Real and imaginary part of the error'
       write (6,*)

       write (6,101) None,error1r,error1i
       write (6,101) Ntwo,error2r,error2i
       write (6,101) Nthr,error3r,error3i

c---
c verify the real eigenvalue
c---

       write (6,*)
       write (6,*) "Verified real eigenvalue"
       write (6,*)

       a11 = a11-x1
       a22 = a22-x1
       a33 = a33-x1
       Det = a11*(a22*a33-a23*a32)
     +      -a21*(a12*a33-a13*a32)
     +      +a31*(a12*a23-a13*a22)
       a11 = a11+x1
       a22 = a22+x1
       a33 = a33+x1
       write (6,101) None,Det

c------------------
c  three real roots
c------------------

      else       ! three real roots

       write (6,*)
       write (6,*) " Three real roots"
       write (6,*)
       write (6,*) 'Roots'
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

c---
c verify the eigenvalues
c---

       write (6,*)
       write (6,*) "Verified eigenvalues"
       write (6,*)

       a11 = a11-x1
       a22 = a22-x1
       a33 = a33-x1
       Det = a11*(a22*a33-a23*a32)
     +      -a21*(a12*a33-a13*a32)
     +      +a31*(a12*a23-a13*a22)
       a11 = a11+x1
       a22 = a22+x1
       a33 = a33+x1
       write (6,101) None,Det
       a11 = a11-x2
       a22 = a22-x2
       a33 = a33-x2
       Det = a11*(a22*a33-a23*a32)
     +      -a21*(a12*a33-a13*a32)
     +      +a31*(a12*a23-a13*a22)
       a11 = a11+x2
       a22 = a22+x2
       a33 = a33+x2
       write (6,101) None,Det
       a11 = a11-x3
       a22 = a22-x3
       a33 = a33-x3
       Det = a11*(a22*a33-a23*a32)
     +      -a21*(a12*a33-a13*a32)
     +      +a31*(a12*a23-a13*a22)
       a11 = a11+x3
       a22 = a22+x3
       a33 = a33+x3
       write (6,101) None,Det

c-----------
      end If
c-----------

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 0 to quit"
      write (6,*) " 1 to repeat"
      write (6,*) " -----------"
      read  (5,*) more
      If(more.eq.1) Go to 98

c-----
c Done
c-----

      write (6,*) " Thank you for your business"

 101  Format (1x,i1,2(1x,f15.10))

      stop
      end
