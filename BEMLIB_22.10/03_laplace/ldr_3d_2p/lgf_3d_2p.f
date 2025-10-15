      subroutine lgf_3d_2p
     +                     (Iopt
     +                     ,Method
     +                     ,x,y,z
     +                     ,x0,y0,z0
     +                     ,a11,a12,a21,a22
     +                     ,b11,b12,b21,b22
     +                     ,ew,area
     +                     ,Max1,Max2,Max3
     +                     ,G
     +                     ,Gx,Gy,Gz
     +                     )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement.
c----------------------------------------

c--------------------------------------
c Computation of the doubly-periodic Green's function
c of Laplace's equation in three dimensions
c
c The point sources are located in a plane
c that is parallel to the xy plane.
c and perpendicular to the z axis
c
c One of the point sources is located at (x0, y0, z0)
c
c The Green's function is computed by two methods
c according to the value of the switch "Method"
c
c SYMBOLS:
c --------
c
c Method = 1 for the straight Fourier series method
c Method = 2 for the fast summation method of
c                Hautman and Klein
c
c x, y, z:	coordinates of the field point
c x0, y0, z0:	coordinates of any one of the point sources
c
c Gx, Gy, Gz:  components of grad(G)
c
c Iopt =  1  compute only G
c      ne 1  compute G and grad(G)
c
c--------------------------------------

      Implicit Double Precision (a-h,k,o-z)

c---
c constants
c---

      pi   = 3.14159 265358
      pi2  = 2.0*pi
      pi4  = 4.0*pi
      srpi = sqrt(pi)

      toe =  3.0D0/8.0D0
      tot =  3.0D0/2.0D0
      foe = 15.0D0/8.0D0

c--------
c prepare
c--------

      dxx = x-x0
      dyy = y-y0
      dzz = z-z0

c-----------
c initialize
c-----------

      G  = 0.0        ! Green's function
      Gx = 0.0        ! and gradient
      Gy = 0.0
      Gz = 0.0

      If(method.eq.1) Go to 83   ! straight Fourier series

c-----------------------------------------
c  Sum in real space of the 3D regularized
c  GF over the 2D lattice
c-----------------------------------------

c-----
c compute the accelerated sum of 1/r
c-----

      fc = 2.0/srpi

      dz  = dzz 
      dzs = dz**2
      dzc = dzs*dz
      dzq = dzs**2

      Do 1 i1 = -Max3,Max3
       Do 1 i2 = -Max3,Max3

        dx = dxx - i1*a11 - i2*a21 
        dy = dyy - i1*a12 - i2*a22

        If(abs(dx)+abs(dy)
     +            +abs(dz).lt.0.00000001) Go to 1    ! field point on
                                                     ! a singularity
        sns = dx**2+dy**2
        sn  = sqrt(sns)

        rns = sns+dzs
        rn  = sqrt(rns)

        en  = dzs/sns
        ens = en**2

        G = G + 1.0/rn-(1.0-0.5*en+toe*ens)/sn

c-------------------------
        ! compute grad(G)

        If(Iopt.ne.1) then

          rnc = rn*rns
          snc = sn*sns
          tmp = (1.0-1.5*en+foe*ens)/snc
          Gx  = Gx-dx/rnc+dx*tmp
          Gy  = Gy-dy/rnc+dy*tmp
          Gz  = Gz-dz/rnc-(-1.0+tot*en)*dz/snc

        End If
c------------------------

  1   Continue

c-----------------------
c compute sums of 1/s**m
c in real space
c-----------------------

c     test0 = 0.0     ! testing invariance with respect to ewald
c     test1 = 0.0
c     test2 = 0.0

      sum0 = 0.0
      sum1 = 0.0
      sum2 = 0.0

      If(Iopt.ne.1) then
       sum0x = 0.0
       sum0y = 0.0
       sum1x = 0.0
       sum1y = 0.0
       sum2x = 0.0
       sum2y = 0.0
      End If

      Do 2 i1 = -Max1,Max1
       Do 2 i2 = -Max1,Max1

        dx = dxx - i1*a11 - i2*a21 
        dy = dyy - i1*a12 - i2*a22

        If((abs(dx)+abs(dy)).lt.0.00000001) Go to 2    ! field point on
                                                       ! a singularity
        sns = dx**2+dy**2    ! square of s_n
        sn  = sqrt(sns)      ! s_n
        snc = sns*sn         ! thrid power of s_n
        snp = sns*snc        ! fifth power of s_n

        w  = ew*sn
        ws = w**2
        wq = ws**2
        we = wq*ws

c---
c compute the error function
c erf(w)
c---

        T = 1.0/(1.0+0.5*w)          ! complementary error function
        ERFCC=    
     +  T*EXP(-W**2-1.26551223+T*(1.00002368+T*(.37409196+
     +  T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+
     +  T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
        erfc = 1.0-erfcc

        expp = exp(-ws)

        h0 = erfc
        h1 = h0-fc*w*expp*(1.0+2.0*ws)
        h2 = h0-fc*w*expp*(9.0+6.0*ws-4.0*wq+8.0*we)/9.0

        sum0 = sum0+(1.0-h0)/sn
        sum1 = sum1+(1.0-h1)/snc
        sum2 = sum2+(1.0-h2)/snp

c-------------------------
        ! compute grad(G)

        If(Iopt.ne.1) then

         snh = sns*snp    ! seventh power of s_n

         wq = ws**2
         dh0dw = fc*expp
         dh1dw = dh0dw *  4.0*ws*(-1.0+ws)
         dh2dw = dh0dw * 16.0*wq*( 2.0-4.0*ws+wq)/9.0

         tmp   = (1.0-h0+w*dh0dw)/snc
         sum0x = sum0x - dx*tmp
         sum0y = sum0y - dy*tmp

         tmp   = (3.0-3.0*h1+w*dh1dw)/snp
         sum1x = sum1x - dx*tmp
         sum1y = sum1y - dy*tmp

         tmp   = (5.0-5.0*h2+w*dh2dw)/snh
         sum2x = sum2x - dx*tmp
         sum2y = sum2y - dy*tmp

        End If
c------------------------- 

 2    Continue

c     write (6,*) " Exited sum in 2D real space"

      G = G + sum0 - 0.5*dzs*sum1  + toe*dzq*sum2

c-------------------------        ! compute grad(G)
      If(Iopt.ne.1) then
        Gx = Gx + sum0x - 0.5*dzs*sum1x + toe*dzq*sum2x
        Gy = Gy + sum0y - 0.5*dzs*sum1y + toe*dzq*sum2y
        Gz = Gz - dz*sum1 + tot*dzc*sum2
      End If
c-------------------------

c---------------------------------------
c     test0 = test0 + sum0   ! testing invariance with respect to ewald
c     test1 = test1 - 0.5*dzs*sum1 
c     test2 = test2 + toe*dzq*sum2
c---------------------------------------

c-------------------------------
c END OF SUMMATION IN REAL SPACE
c-------------------------------

 83   Continue

c-------------------------
c SUM IN WAVENUMBER SPACE
c-------------------------

c-------------------
c     p0 = 0.0   ! testing invariance with respect to ew
c     p1 = 0.0
c     p2 = 0.0
c     testt0 = 0.0   
c     testt1 = 0.0
c     testt2 = 0.0
c-------------------

      p  = 0.0
      px = 0.0
      py = 0.0
      pz = 0.0

      Do 3 i1=-Max2,Max2
       Do 3 i2=-Max2,Max2

        k1 = i1*b11 + i2*b21 
        k2 = i1*b12 + i2*b22
 
        ks = k1**2 + k2**2
        k  = sqrt(ks)

        If(k.le.0.000001) Go to 3   ! skip the zero point

c-------------------------
      If(Method.eq.1) then          ! straight Fourier series
c-------------------------

        rho = abs(dzz)*k

        A = exp(-rho)/k

        If(Iopt.ne.1) then     ! for grad(G)
         B =  k1*A
         C =  k2*A
         D = -k *A
        End If

c-------------------------
      Else If(Method.eq.2) then     ! accelerated
c-------------------------

        w = k/(2.0*ew)
        T = 1.0/(1.0+0.5*w)          ! complementary error function
        ERFCC=    
     +  T*EXP(-W**2-1.26551223+T*(1.00002368+T*(.37409196+
     +  T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+
     +  T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
 
        fq0 =  erfcc/k 
        fq1 = -ks*fq0
        fq2 = -ks*fq1/9.0

        A = fq0 - 0.5*dzs*fq1 + toe*dzq*fq2

        If(Iopt.ne.1) then    ! for grad(G)
         B = k1*A
         C = k2*A
         D =  - dz*fq1 + tot*dzc*fq2
        End If

c-----------------------------------
c       testt0 = testt0 + fq0
c       testt1 = testt1 - 0.5*dzs*fq1 
c       testt2 = testt2 + toe*dzq*fq2
c-----------------------------------

c-----------
      End If
c-----------

       arg = k1*dxx+k2*dyy
       cs  = cos(arg)

       p = p + A*cs

       If(Iopt.ne.1) then         ! for grad(G)
        sn = sin(arg)
        px = px - B*sn
        py = py - C*sn
        pz = pz + D*cs
       End If

c---------------------
c       p0 = p0+fq0*cs
c       p1 = p1+fq1*cs
c       p2 = p2+fq2*cs
c---------------------

c---
c end of computing grad(G)
c---

  3   Continue

c-------------------------------
c END OF SUM IN RECIPROCAL SPACE
c-------------------------------

c------------------------
c straight Fourier series
c------------------------

      If(method.eq.1) then

        fc = 0.5/area

        G = fc*p
        G = G - fc*abs(dzz)      ! regularize by adding a linear function

        If(Iopt.ne.1) then       ! for grad(G)

         Gx = fc*px              ! normalize
         Gy = fc*py

         If(dzz.lt.0) pz = -pz
         Gz = fc*pz - fc*dzz/abs(dzz)               ! regularize 

        End If

      End If

c-----------------
c expedited method
c-----------------

      If(method.eq.2) then

       prf = pi2/area

       p = p - 1.0/(ew*srpi)  ! contribution from zero wave number
       G = G + prf*p          ! add the sums in real and reciprocal space

       G = G/pi4              ! normalize

       If(Iopt.ne.1) then     ! for grad(G)

        Gx = Gx + prf*px      ! add the sums in real and reciprocal space
        Gy = Gy + prf*py
        Gz = Gz + prf*pz

        Gx = Gx/pi4           ! normalize
        Gy = Gy/pi4
        Gz = Gz/pi4

       End If

c-------------------------------
c      testt0 = testt0  - 2.0*srpi/ew
c      p0 = p0 - 2.0*srpi/ew
c      p0 = fc*p0
c      p1 = fc*p1
c      p2 = fc*p2
c      write (6,*) sum0,sum1,sum2
c      write (6,*) p0,p1,p2
c      write (6,*) p0+sum0,p1+sum1,p2+sum2
c
c  test0  = test0 + testt0*fc      ! should be independent of xi
c  test1  = test1 + testt1*fc
c  test2  = test2 + testt2*fc
c  write (6,*) test0,test1,test2
c-------------------------------

      End If

c---
c Done
c---

  100 Format (9(1x,f10.6))

      Return
      End

c================================================

      subroutine ewald_3d_2p
     +                      (a11,a12,a21,a22
     +                      ,b11,b12,b21,b22
     +                      ,ew,area
     +                      )

c---------------------------------------------
c Computes the reciprocal lattice base vectors
c and the optimal value of the parameter xi
c---------------------------------------------

      Implicit Double Precision (a-h,o-z)

c---
c constants
c---

      pi   = 3.14159 265358
      pi2  = 2.0*pi
      srpi = sqrt(pi)

c---
c unit cell area
c---

      area = a11*a22-a21*a12

c-----------------------------------------
c lattice base vectors in wave number space
c-----------------------------------------

      f   = pi2/area

      b11 =  f*a22
      b12 = -f*a21
      b21 = -f*a12
      b22 =  f*a11

c     write (6,*) " Reciprocal Lattice "
c     write (6,*) " ------------------ "
c     write (6,104) b11,b12
c     write (6,104) b21,b22

      ew = 0.5*srpi/sqrt(area)

c     write (6,114) area,ew

c---
c Done
c---

 104  Format (6(1x,f10.5))
 114  Format (2x," area=",f10.5," ewald parameter =",f10.5)

      Return
      End
