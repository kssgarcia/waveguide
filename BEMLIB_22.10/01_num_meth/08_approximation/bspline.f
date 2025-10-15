      program bspline

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c------------------------------------------------
c    This program accompanies the book:
c
c           C. Pozrikidis
c Numerical Computation in Science and Engineering
c      Oxford University Press
c------------------------------------------------

c----------------------------------------
c B-spline approximation of kth degree
c
c Section 8.4
c
c SYMBOLS:
c -------
c
c np:     number of intervals
c
c npl:    number of divisions for plotting
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension x(-200:200),f(-200:200)
      Dimension A(0:20,-100:100)

      Integer p,q

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi

      null = 0

c---
c prepare
c---

      write (6,*)
      write (6,*) " Please enter the number of intervals"
      write (6,*) " ------------------------------------"
      read  (5,*) n

      write (6,*)
      write (6,*) " Please enter the number of plotting points"
      write (6,*) " between the influence points"
      write (6,*) " ----------------------------"
      read  (5,*) npl

c--------
c prepare
c--------

      nprint = n*npl+1      ! total number of printed points

      open (2,file="PLOTDAT")

c-------------------
c Generate the knots
c-------------------

      Do i=-n,2*n
        tmp  = (i-1.0D0)/n
        arg  = pi2*tmp              ! example
        x(i) = tmp-0.1D0*Dsin(arg)    ! example
        arg  = pi2*x(i)             ! example
        f(i) = Dsin(arg)             ! example
      End Do

      write (6,*) 
      write (6,*) " knots:"
      write (6,*) 

      write (2,100) n+1

      Do i=1,n+1
        write (2,100) i,x(i),f(i)
        write (6,100) i,x(i),f(i)
      End Do

c---

  98  Continue

      write (6,*) 
      write (6,*) " Please choose the degree of the spline: k"
      write (6,*) 
      write (6,*) " Enter 0 to quit"
      write (6,*) "       1 for linear "
      write (6,*) "       2 for quadratic "
      write (6,*) "       etc "
      write (6,*) " ----------------------"
      read  (5,*) k

      If(k.eq.0) Go to 99

c--------------------------------
c Will produce B splines
c according to equations (8.4.19) and (8.4.20)
c--------------------------------

      If(mod(k,2).ne.0) then
        write (6,*) "Confirmed that k is Odd"
        kh  = (k+1)/2
        ieo = 1
      Else
        write (6,*) "Confirmed that k is Even"
        kh  = k/2
        ieo = 2
      End If

c---------------------------------------
c compute approximation at nprint points
c---------------------------------------

      icount = 1

      write (6,*) 
      write (6,*) " Approximation and exact value"
      write (6,*) 
      write (2,100) nprint

      Do 1 i=1,n      ! run over the knots

      If(ieo.eq.1) then    ! odd-degree splines
        imin = i-kh+1
        imax = i+kh
      Else                 ! even-degree splines
        imin = i-kh
        imax = i+kh
      End If

c-----------
c initialize
c-----------

      Do p=0,k
        Do j=-k,n+k
          A(p,j) = 0.0D0
        End Do
      End Do

      A(0,i) = 1.0D0

      ippmax = npl
      If(i.eq.n) ippmax = npl+1    ! to get the very last point

      Dx = (x(i+1)-x(i))/npl
      xx = x(i)

c---
c Cox-deBoor recursion
c---

      Do ipp=1,ippmax    ! loop over plotting points

        Do p=1,k
         Do q=i-p,i
         A(p,q) = (xx-x(q))    *A(p-1,q)  /(x(p+q)  -x(q))
     +          + (x(q+p+1)-xx)*A(p-1,q+1)/(x(p+1+q)-x(q+1))
         End Do
        End Do

        ff = 0.0               ! approximation
        If(ieo.eq.1) Then
         Do j=imin,imax
           ff = ff + f(j)*A(k,j-kh)
         End Do
        Else
         Do j=imin,imax
           ff = ff + 0.5D0*(f(j)+f(j+1))*A(k,j-kh)
         End Do
        End If

        exact = Dsin(xx*pi2)

        write (6,100) icount,xx,ff,exact
        write (2,100) icount,xx,ff,exact

        icount = icount+1
        xx = xx+Dx

      End Do                    ! end of loop over points

  1   Continue

c-------
c repeat
c-------

      Go to 98

c-----
c Done
c-----

  99  Continue

      write (2,100) null
      write (6,100) null

c-----
c Rest
c-----

 100  Format (1x,i3,5(1x,f10.5))

      Stop
      End
