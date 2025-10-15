      subroutine elm_line
     +
     +  (N
     +  ,ratio
     +  ,X1,Y1
     +  ,X2,Y2
     +  ,sinit
     +  ,Isym
     +  ,Xe,Ye,se
     +  ,Xm,Ym,sm
     +  )

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-----------------------------------------------------
c Disretization of a line segment into N elements
c
c  X1,Y1: coordinates of the first point
c  X2,Y2: coordinates of the last point
c
c  ratio: 
c
c     If Isym = 0, ratio of length of LAST to FIRST element
c     If Isym = 1, ratio of length of MID  to FIRST element
c
c  alpha: geometric factor ratio
c
c  sinit: specified arc length at (X1, Y1)
c
c  se: arc length at the element end-nodes
c  sm: arc length at the element mid-nodes
c
c  xe,ye: end nodes
c  xm,ym: mid nodes
c
c-----------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension Xe(100),Ye(100),se(100)
      Dimension Xm(100),Ym(100),sm(100)

c------------
c one element
c------------

      If(N.eq.1) then

       xe(1) = X1 
       ye(1) = Y1 

       xe(2) = X2 
       ye(2) = Y2 

       se(1) = sinit
       se(2) = se(1)+Dsqrt( (X2-X1)**2+(Y2-Y1)**2)

       Go to 99

      End If

c--------
c prepare
c--------

      N1 = N+1

c--------------------
c biased distribution
c--------------------

      If(Isym.eq.0) then

      If(ratio.eq.1.000) then
        alpha  = 1.0D0
        factor = 1.0D0/N
      Else
        texp   = 1.0D0/(N-1.0D0)
        alpha  = ratio**texp
        factor = (1.0D0-alpha)/(1.0D0-alpha**N)
      End If

      deltax = (x2-x1) * factor   ! x length of first element
      deltay = (y2-y1) * factor   ! y length of first element 

      Xe(1) = X1     ! first point
      Ye(1) = Y1
      se(1) = sinit

      Do i=2,N1
        Xe(i)  = Xe(i-1)+deltax
        Ye(i)  = Ye(i-1)+deltaY
        se(i)  = se(i-1)+Dsqrt(deltax**2+deltay**2)
        deltax = deltax*alpha
        deltaY = deltay*alpha
      End Do

      Go to 99

      End If

c----------------------
c symmetric distribution
c even number of points
c----------------------

      If(mod(N,2).eq.0) then
   
      Xh = 0.5D0*(X1+X2)      ! mid-point
      Yh = 0.5D0*(Y1+Y2)

      Nh  = N/2
      Nh1 = Nh+1

      If(ratio.eq.1.000) then
        alpha  = 1.0D0
        factor = 1.0D0/Nh
      Else
        texp   = 1.0D0/(Nh-1.0)
        alpha  = ratio**texp
        factor = (1.0D0-alpha)/(1.0D0-alpha**Nh)
      End If

      deltax = (xh-x1) * factor   ! x length of first element
      deltay = (yh-y1) * factor   ! y length of first element 

      Xe(1) = X1    ! first point
      Ye(1) = Y1
      se(1) = sinit

      Do i=2,Nh1
        Xe(i)  = Xe(i-1)+deltax
        Ye(i)  = Ye(i-1)+deltaY
        se(i)  = se(i-1)+sqrt(deltax**2+deltay**2)
        deltax = deltax*alpha
        deltaY = deltay*alpha
      End Do

      deltax = deltax/alpha
      deltaY = deltay/alpha

      Do i=Nh1+1,N1
        Xe(i)  = Xe(i-1)+deltax
        Ye(i)  = Ye(i-1)+deltaY
        se(i)  = se(i-1)+sqrt(deltax**2+deltay**2)
        deltax = deltax/alpha
        deltay = deltay/alpha
      End Do

      Go to 99

      End If

c-----------------------
c symmetric distribution
c odd number of points
c-----------------------

      If(ratio.eq.1.000) then
        alpha  = 1.0D0
        factor = 1.0D0/N1
      Else
        texp   = 2.0D0/(N-1.0D0)
        alpha  = ratio**texp
        tmp1   = 0.50D0*(N+1.0D0)
        tmp2   = 0.50D0*(N-1.0D0)
        factor = (1.0D0-alpha)/(2.0D0-alpha**tmp1-alpha**tmp2)
      End If

      deltax = (x2-x1) * factor   ! x length of first element
      deltay = (y2-y1) * factor   ! y length of first element 

      Xe(1) = X1     ! first point
      Ye(1) = Y1
      se(1) = sinit

      Do i=2,(N+3)/2
        Xe(i)  = Xe(i-1)+deltax
        Ye(i)  = Ye(i-1)+deltaY
        se(i)  = se(i-1)+sqrt(deltax**2+deltay**2)
        deltax = deltax*alpha
        deltaY = deltay*alpha
      End Do

      deltax = deltax/(alpha**2)
      deltay = deltay/(alpha**2)

      Do i=(N+5)/2,N1
        Xe(i)  = Xe(i-1)+deltax
        Ye(i)  = Ye(i-1)+deltaY
        se(i)  = se(i-1)+sqrt(deltax**2+deltay**2)
        deltax = deltax/alpha
        deltaY = deltay/alpha
      End Do

c-----
c Done
c-----

  99  Continue

c-------------------
c compute mid-points
c-------------------

      Do i=1,N
       xm(i) = 0.5D0*(xe(i)+xe(i+1))
       ym(i) = 0.5D0*(ye(i)+ye(i+1))
       sm(i) = 0.5D0*(se(i)+se(i+1))
      End Do

c-----
c Done
c-----

      write (6,*) "elm_line: Geometric ratio: ",alpha

      Return
      End
