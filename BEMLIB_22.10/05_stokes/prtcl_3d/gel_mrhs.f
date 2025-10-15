      subroutine gel_mrhs 
     +
     +   (n
     +   ,a
     +   ,nrhs,rhs
     +   ,x
     +   ,Isym,Iwlpvt
c    +   ,l,u
     +   ,det
     +   ,Istop
     +   )

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-------------------------------------------------
c  Gauss elimination with nrhs right-hand sides,
c  with option of row pivoting.
c
c  This subroutine returns solution vectors, upper
c  and lower triangular factors, flag for completion,
c  and the determinant.
c---------------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c              C. Pozrikidis
c Numerical Computation in Science and Engineering
c         Oxford University Press
c------------------------------------------------

c________________   external variables   _______________
c
c  a ...... square matrix
c  n ...... size (rows/columns) of matrix a
c  nrhs.... number of right-hand sides
c  rhs .... columns are right hand side vectors (e.g. b, as in Ax=b)
c  x ...... columns are solution vectors
c  Isym ... flag for symmetry of matrix a (1 = symmetric)
c  Iwlpvt.. 0 for no pivoting, 1 for pivoting
c  eps..... tolerance to identify a singular matrix
c  tol..... tolerance for the residuals
c  l ...... lower triangular matrix
c  u ...... upper triangular matrix
c  det .... determinant (det(a)=det(l)*det(u))
c   Istop... flag: Istop = 1 if something is wrong
c
c________________   internal variables   _______________
c       
c   pivot .. absolute value of pivot candidates
c   ipv..... location of pivotal element
c   icount . counter for number of row interchanges
c_______________________________________________________
c
c---------------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(1600,1600),rhs(1600,10),c(1600,1610),u(1600,1600)
      Double Precision l(1600,1600)
      Dimension x(1600,10)

      Parameter (eps=0.000001,tol=0.0000001)

c----------
c flags etc
c----------

      Istop =0

      Icount=0

c--------
c prepare
c--------

      na = n-1
      n1 = n+1
      ntot = n+nrhs

c-------------------
c initialize l and c
c-------------------

      Do i=1,n
       Do j=1,n
         l(i,j) = 0.0D0
         c(i,j) = a(i,j)
       End Do
       Do j=1,nrhs
         c(i,n+j) = rhs(i,j)
       End Do
      End Do

c---------------------
c begin row reductions
c---------------------

      Do 1 m=1,na            ! outer loop for working row

      ma = m-1
      m1 = m+1

      If(Iwlpvt.ne.1) Go to 97    ! skip pivoting

c---
c Pivoting
c begin by searching column i for largest element
c---

      ipv  = m
      pivot = Dabs(c(m,m))

      Do j=m1,n
        If(Dabs(c(j,m)).gt.pivot) then
         ipv   = j
         pivot = Dabs(c(j,m))
        End If
      End Do

      If(pivot.lt.eps) then
        write (6,*) 
        write (6,*) " gauss_mrhs: trouble at station 1"
        write (6,*) 
        Istop = 1
        Return
      End If

c---
c  switch the working row with the row containing the 
c  pivot element (also switch rows in l)
c---

      If(ipv.ne.m) then
        Do j=m,ntot
         save    = c(m,j)
         c(m,j)  = c(ipv,j)
         c(ipv,j)= save
        End Do

        Do j=1,ma
         save     = l(m,j)
         l(m,j)   = l(ipv,j)
         l(ipv,j) = save
        End Do

        icount=icount+1 

      End If

 97   Continue

c---
c reduce column i beneath element c(m,m)
c---

      Do 2 i=m1,n

c-------------------------
        If(Isym.eq.1) then
c-------------------------
          l(i,m) = c(m,i)/c(m,m)
          Do j = i,ntot
           c(i,j)=c(i,j)-l(i,m)*c(m,j)
          End Do
c-----------
        Else
c-----------
          l(i,m) = c(i,m)/c(m,m)
          c(i,m) = 0.0
          Do j=m1,ntot
           c(i,j)=c(i,j)-l(i,m)*c(m,j)
          End Do
c-------------
        End If
c-------------

 2    Continue

c---
c fill in the leading zeros in the matrix c
c (not necessary)
c---

c      Do j=1,m
c        c(m+1,j)=0.0D0
c      End Do

 1    Continue         ! end of outer loop for working row

c---
c check the last diagonal element for singularity
c---

      If(abs(c(n,n)).lt.eps) then
        write (6,*) " Trouble at station 2 of Gauss elimin."
        Istop = 1
        Return
      End If

c----------------------
c complete the matrix l
c----------------------

      Do i=1,n
        l(i,i)=1.0D0
      End Do

c--------------------
c define the matrix u
c--------------------

      Do i =1,n
        Do j =1,n
         u(i,j) = c(i,j)
        End Do
      End Do

c---
c perform back-substitution to solve the reduced system
c using the upper triangular matrix c
c---

      Do ll=1,nrhs

        x(n,ll) = c(n,n+ll)/c(n,n)

        Do i=n-1,1,-1
          sum=c(i,n+ll)
          Do j=i+1,n
          sum=sum-c(i,j)*x(j,ll)
          End Do
          x(i,ll) = sum/c(i,i)
        End Do

      End Do

c---
c compute determinant as det(a)= +- det(l)*det(u)
c---
      Det=1.0D0

      Do i=1,n
        det=det*c(i,i)
      End Do

      If(Iwlpvt.eq.1) then
        write (6,*)
        write (6,*) " Number of row interchanges : ",Icount
        write (6,*)
        Do i=1,icount
          det = -det
        End Do
      End If

c---
c compute residuals
c---

      Do ll=1,nrhs
       Do i =1,n
        sum = rhs(i,ll)
        Do j=1,n
        sum = sum - a(i,j)*x(j,ll)
        End Do
        If(Dabs(sum).gt.tol) then
          Istop = 1
          write (6,100) i,sum
        End If
       End Do
      End Do

c-----
c done
c-----

  100 Format (1x,i4,f15.10)
  101 Format(16(16(1x,f5.3),/))

      Return
      End
