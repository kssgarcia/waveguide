      subroutine elm_arc
     +
     +  (N
     +  ,ratio
     +  ,Xcnt,Ycnt
     +  ,Radius
     +  ,angle1,angle2
     +  ,sinit
     +  ,Isym
     +  ,Xe,Ye,Te,se
     +  ,Xm,Ym,Tm,sm
     +  )

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c----------------------------------------------------------------
c  Disretization of a circular arc into N elements
c
c  SYMBOLS:
c  -------
c
c  Xcnt, Ycnt:	center of the arc
c  ratio:	radius of the arc
c
c  ratio: If Isym = 0, ratio of length of LAST to FIRST element
c         If Isym = 1, ratio of length of MID  to FIRST element
c
c  Xe, Ye, Te:	coordinates and polar angle for the end-nodes
c  Xm, Ym, Tm:	coordinates and polar angle for the mid-nodes
c
c  sinit: assigned arc length at the first point
c----------------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension Xe(100),Ye(100),Te(100),se(100)
      Dimension Xm(100),Ym(100),Tm(100),sm(100)

c------------
c one element
c------------

      if(N.eq.1) then

       xe(1) = X1
       ye(1) = Y1
       te(1) = angle1 
       se(1) = sinit

       xe(2) = X2
       ye(2) = Y2
       te(2) = angle12
       se(2) = se(1)+Dsqrt( (X2-X1)**2+(Y2-Y1)**2)

       Go to 99

      end if

c--------
c prepare
c--------

      N1 = N+1

c---------------------------
c non-symmetric distribution
c---------------------------

      if(Isym.eq.0) then

      texp  = 1.0D0/(N-1.0D0)
      alpha = ratio**texp
 
      if(Dabs(alpha-1.0D0).gt.0.0000001) then
       factor = (1.0D0-alpha)/(1.0D0-alpha**N)
      else
       factor = 1.0D0/(N-1.0D0+1.0D0)
      end if

      deltat = (angle2-angle1) * factor   ! aperture of first element

      te(1) = angle1                      ! first point
      Xe(1) = Xcnt + radius*Dcos(te(1))
      Ye(1) = Ycnt + radius*Dsin(te(1))
      se(1) = sinit

      Do i=2,N1
        Te(i) = Te(i-1)+deltat
        Xe(i) = Xcnt + radius*Dcos(Te(i))
        Ye(i) = Ycnt + radius*Dsin(Te(i))
        se(i) = se(i-1)+radius*Dabs(deltat)
        deltat = deltat*alpha
      End Do

      Go to 99

      end if

c-----------------------
c symmetric distribution
c----------------------

c---
c even number of elements
c---

      if(mod(N,2).eq.0) then

      angleh = 0.5D0*(angle1+angle2)      ! mid-point

      Nh  = N/2
      Nh1 = Nh+1

      texp  = 1.0D0/(Nh-1.0D0)
      alpha = ratio**texp

      if(abs(alpha-1.0D0).gt.0.0000001) then
       factor = (1.0D0-alpha)/(1.0D0-alpha**Nh)
      else
       factor = 1.0D0/(Nh-1.0D0+1.0D0)
      end if

      deltat = (angleh-angle1) * factor   ! aperture of first element

      te(1) = angle1                      ! first point
      Xe(1) = Xcnt + radius*Dcos(te(1))
      Ye(1) = Ycnt + radius*Dsin(te(1))
      se(1) = sinit

      Do i=2,Nh1                          ! up to mid-point
        te(i) = te(i-1)+deltat
        Xe(i) = Xcnt + radius*Dcos(te(i))
        Ye(i) = Ycnt + radius*Dsin(te(i))
        se(i) = se(i-1)+radius*Dabs(deltat)
        deltat = deltat*alpha
      End Do

      deltat = deltat/alpha

      Do i=Nh1+1,N1                        ! reflect
        te(i) = te(i-1)+deltat
        Xe(i) = Xcnt + radius*Dcos(te(i))
        Ye(i) = Ycnt + radius*Dsin(te(i))
        se(i) = se(i-1)+radius*Dabs(deltat)
        deltat = deltat/alpha
      End Do

      Go to 99

      end if

c---
c odd number of elements
c---

      texp  = 2.0D0/(N-1.0D0)
      alpha = ratio**texp

      if(Dabs(alpha-1.0D0).gt.0.0000001) then
       tmp1   = (N+1.0)/2.0D0
       tmp2   = (N-1.0)/2.0D0
       factor = (1.0D0-alpha)/(2.0D0-alpha**tmp1-alpha**tmp2)
      else
       factor = 1.0/N
      end if

      deltat = (angle2-angle1) * factor   ! aperture of first element

      te(1) = angle1                      ! first point
      Xe(1) = Xcnt + radius*Dcos(te(1))
      Ye(1) = Ycnt + radius*Dsin(te(1))
      se(1) = sinit

      Do i=2,(N+3)/2                      ! up to mid point + 1
        Te(i) = Te(i-1)+deltat
        Xe(i) = Xcnt + radius*Dcos(te(i))
        Ye(i) = Ycnt + radius*Dsin(te(i))
        se(i) = se(i-1)+radius*Dabs(deltat)
        deltat = deltat*alpha
      End Do

      deltat = deltat/(alpha**2)

      Do i=(N+5)/2,N1
        Te(i) = Te(i-1)+deltat
        Xe(i) = Xcnt + radius*Dcos(te(i))
        Ye(i) = Ycnt + radius*Dsin(te(i))
        se(i) = se(i-1)+radius*Dabs(deltat)
        deltat = deltat/alpha
      End Do

c-----------------------------------

  99  Continue

c-----------------------
c compute the mid-points
c-----------------------

      Do i=1,N
       tm(i) = 0.5D0*(te(i)+te(i+1))
       sm(i) = 0.5D0*(se(i)+se(i+1))
       Xm(i) = Xcnt + radius*Dcos(tm(i))
       Ym(i) = Ycnt + radius*Dsin(tm(i))
      End Do

c-----
c Done
c-----

      Return
      End
