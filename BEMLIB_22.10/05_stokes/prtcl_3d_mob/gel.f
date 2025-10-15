      subroutine gel
     +
     +   (n
     +   ,a,rhs,x
     +   ,Isym,Iwlpvt
c    +   ,l,u
     +   ,det
     +   ,Istop
     +   )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c C. Pozrikidis
c
c Numerical Computation in Science and Engineering
c
c Oxford University Press
c
c 1998
c------------------------------------------------

c-----------------------------------------------
c  Gauss elimination with option of row pivoting. 
c
c  Algorithm 3.5.1
c
c  This subroutine returns:
c
c    the solution vector,
c    the and lower triangular matrix factors L and U,
c    a flag for completion
c    the determinant
c
c  SYMBOLS:
c  --------
c
c   a ...... square matrix
c   n ...... size (rows/columns) of matrix a
c   rhs .... right hand side vector (e.g. b, as in Ax=b)
c   c....... extended matrix
c   x ...... solution vector
c
c   Isym ... flag = 1 if a is symmetric; 0 if nonsymmetric
c   Iwlpvt.. 0 for no pivoting, 1 for pivoting
c
c   eps..... tolerance to identify a singular matrix
c   tol..... tolerance for the residuals
c
c   l ...... lower triangular matrix
c   u ...... upper triangular matrix
c   det .... determinant (det(a) = +- det(l)*det(u))
c
c   Istop... flag: Istop = 1 if something is wrong
c
c   pivot .. absolute value of pivot candidates
c   ipv..... location of pivotal element
c   Icount.. counter for number of row interchanges
c
c-------------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension A(1542,1542),rhs(1542),c(1542,1543),u(1542,1542)
      Double Precision l(1542,1542)
      Dimension x(1542)

      Parameter (eps=0.000001,tol=0.0000001)

c-----------
c initialize
c-----------

      Istop  = 0
      Icount = 0     ! counts row interchanges

c--------
c prepare
c--------

      na = n-1
      n1 = n+1

c-------------------
c Initialize l and c
c-------------------

      Do i=1,n
       Do j=1,n
         l(i,j) = 0.0D0
         c(i,j) = a(i,j)
        End Do
        c(i,n1) = rhs(i)
      End Do

c---------------------
c Begin row reductions
c---------------------

      Do 1 m=1,na           ! outer loop for working row

      ma = m-1
      m1 = m+1
        
c-----------------------------
c Pivoting module
c
c begin by searching column i 
c for largest element
c----------------------------

      If(Iwlpvt.ne.1) Go to 97   ! skip pivoting module
   
      ipv  = m
      pivot = abs(c(m,m))
 
      Do j=m1,n
       If(abs(c(j,m)).gt.pivot) then
        ipv   = j
        pivot = abs(c(j,m))
       End If
      End Do

      If(pivot.lt.eps) then
        write (6,*)
        write (6,*) " Trouble at station 1 of Gauss Elimination"
        write (6,*)
        Istop = 1
        Return
      End If

c--------------------------------------
c switch the working row with
c the row containing the pivot element
c
c also switch rows in l
c--------------------------------------
 
      If(ipv.ne.m) then

       Do j=m,n1
         save     = c(m,j)
         c(m,j)   = c(ipv,j)
         c(ipv,j) = save
       End Do

       Do j=1,ma
         save     = l(m,j)
         l(m,j)   = l(ipv,j)
         l(ipv,j) = save
       End Do

       Icount = Icount+1 

      End If

 97   Continue        ! End of pivoting module

c---------------------------------------
c reduce column i beneath element c(m,m)
c---------------------------------------

      Do 2 i=m1,n

       If(Isym.eq.1) then        ! symmetric matrix

         l(i,m) = c(m,i)/c(m,m)
         c(i,m) = 0.0D0

         Do j=i,n1
          c(i,j) = c(i,j)-l(i,m)*c(m,j)
         End Do

        Else                     ! non-symmetric matrix

         l(i,m) = c(i,m)/c(m,m)
         c(i,m) = 0.0D0

         Do j=m1,n1
          c(i,j) = c(i,j)-l(i,m)*c(m,j)
         End Do

       End If

 2    Continue 

 1    Continue                ! end of outer loop for working row 

c---------------------------------
c check the last diagonal element
c for singularity
c--------------------------------

      If(abs(c(n,n)).lt.eps) then

        write (6,*)
        write (6,*) " gauss: trouble at station 2"
        write (6,*) 
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

      Do i=1,n
        Do j=1,n
         u(i,j) = c(i,j)
        End Do
      End Do

c-----------------------------------
c perform back-substitution to solve
c the reduced system
c using the upper triangular matrix c
c------------------------------------

      x(n) = c(n,n1)/c(n,n)

      Do i=na,1,-1
        sum=c(i,n1)
        Do j=i+1,n
         sum=sum-c(i,j)*x(j)
        End Do
        x(i)=sum/c(i,i)
      End Do

c-----------------------
c compute the determinant as
c
c det(a) = (+-) det(l)*det(u)
c
c-----------------------

      det = 1.0D0

      Do i=1,n
       det=det*c(i,i)
      End Do

      If(Iwlpvt.eq.1) then

        write (6,*) " gel: number of row interchanges : ",Icount

        Do i=1,Icount
         det = -det
        End Do

      End If

c----------------------
c compute the residuals
c----------------------

      Do i=1,n

        sum = rhs(i)
        Do j=1,n
        sum = sum - a(i,j)*x(j)
        End Do

        If(abs(sum).gt.tol) then
          Istop = 1
          write (6,*) " gel: problem in solving the linear system"
          write (6,100) i,sum
        End If

      End Do

c--------
c wrap up
c--------

  100 Format (1x,i4,f15.10)
  101 Format(16(16(1x,f5.3),/))

      Return
      End
