      subroutine lgf_3d_2p
     +
     +  (Iopt
     +  ,method
     +  ,x,y,z
     +  ,x0,y0,z0
     +  ,a11,a12,a21,a22
     +  ,b11,b12,b21,b22
     +  ,ew,area
     +  ,Max1,Max2,Max3
     +  ,G
     +  ,Gx,Gy,Gz
     +  )

c================================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c================================================

c--------------------------------------
c Doubly periodic Green's function
c of Laplace's equation in three dimensions
c
c The point sources are located in a plane
c that is perpendicular to the z axis
c
c One of the point sources is located at (x0, y0, z0)
c
c The Green's function is computed by two methods
c according to the value of the switch "method"
c
c SYMBOLS:
c --------
c
c method = 1: Fourier series method
c method = 2: fast summation method of Hautman & Klein
c
c x, y, z:	coordinates of the field point
c x0, y0, z0:	coordinates of any one of the point sources
c
c G: Green's function
c
c Gx, Gy, Gz:  components of grad(G)
c
c Iopt =  1  compute only G
c      ne 1  compute G and grad(G)
c--------------------------------------

      Implicit Double Precision (a-h,k,o-z)

c----------
c constants
c----------

      pi   = 3.14159 265358 D0
      pi2  = 2.0D0*pi
      pi4  = 4.0D0*pi
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

      G  = 0.0D0        ! Green's function
      Gx = 0.0D0        ! and gradient
      Gy = 0.0D0
      Gz = 0.0D0

c-----------------
c expedited method
c-----------------

      if(method.eq.2) then

c-----------------------------------------
c sum in real space of the 3D regularized
c Green's function over the 2D lattice
c-----------------------------------------

c-----
c  compute the accelerated sum of 1/r
c-----

      fc = 2.0D0/srpi

      dz  = dzz 
      dzs = dz*dz
      dzc = dzs*dz
      dzq = dzs*dzs

      Do i1=-Max3,Max3
       Do i2=-Max3,Max3

        dx = dxx - i1*a11 - i2*a21 
        dy = dyy - i1*a12 - i2*a22

        if((abs(dx)+abs(dy)                      ! field point on
     +            +abs(dz)).gt.0.00000001) then  ! a singularity
                                                
        sns = dx*dx+dy*dy
        sn  = Dsqrt(sns)

        rns = sns+dzs
        rn  = Dsqrt(rns)

        en  = dzs/sns
        ens = en*en

        G = G + 1.0D0/rn - (1.0D0-0.5D0*en+toe*ens)/sn

        ! compute grad(G)

        if(Iopt.ne.1) then
          rnc = rn*rns
          snc = sn*sns
          tmp = (1.0D0-1.5D0*en+foe*ens)/snc
          Gx  = Gx-dx/rnc+dx*tmp
          Gy  = Gy-dy/rnc+dy*tmp
          Gz  = Gz-dz/rnc-(-1.0+tot*en)*dz/snc

        end if
c------------------------

       end if

       End Do
      End Do

c-------------------------------------
c compute sums of 1/s**m in real space
c-------------------------------------

c     test0 = 0.0     ! testing invariance with respect to ewald
c     test1 = 0.0
c     test2 = 0.0

      sum0 = 0.0D0
      sum1 = 0.0D0
      sum2 = 0.0D0

      if(Iopt.ne.1) then
       sum0x = 0.0D0
       sum0y = 0.0D0
       sum1x = 0.0D0
       sum1y = 0.0D0
       sum2x = 0.0D0
       sum2y = 0.0D0
      end if

      Do i1=-Max1,Max1
       Do i2=-Max1,Max1

        dx = dxx - i1*a11 - i2*a21 
        dy = dyy - i1*a12 - i2*a22

        if((abs(dx)+abs(dy)).gt.0.00000001) then    ! field point on
                                                    ! a singularity
        sns = dx**2+dy**2    ! square of s_n
        sn  = Dsqrt(sns)      ! s_n
        snc = sns*sn         ! thrid power of s_n
        snp = sns*snc        ! fifth power of s_n

        w  = ew*sn
        ws = w*w
        wq = ws*ws
        we = wq*ws

c---
c compute the error function
c erf(w)
c---

        T = 1.0/(1.0+0.5*w) 
        ERFCC=             ! complementary error function
     +  T*EXP(-W*W-1.26551223+T*(1.00002368+T*(.37409196+
     +  T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+
     +  T*(1.48851587+T*(-.82215223+T*.17087277)))))))))

        erfc = 1.0D0-erfcc

        expp = exp(-ws)

        h0 = erfc
        h1 = h0- fc * w * expp*(1.0+2.0*ws)
        h2 = h0- fc * w * expp*(9.0+6.0*ws-4.0*wq+8.0*we)/9.0D0

        sum0 = sum0+(1.0D0-h0)/sn
        sum1 = sum1+(1.0D0-h1)/snc
        sum2 = sum2+(1.0D0-h2)/snp

c----------------
c compute grad(G)
c----------------

        if(Iopt.ne.1) then

         snh = sns*snp    ! seventh power of s_n

         wq = ws*ws
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

        end if
c------------------------- 
        end if

       End Do
      End Do
c------------------------- 

c     write (6,*) " Exited sum in 2D real space"

      G = G + sum0 - 0.5D0*dzs*sum1  + toe*dzq*sum2

c-------------------------        ! compute grad(G)
      if(Iopt.ne.1) then
        Gx = Gx + sum0x - 0.5D0*dzs*sum1x + toe*dzq*sum2x
        Gy = Gy + sum0y - 0.5D0*dzs*sum1y + toe*dzq*sum2y
        Gz = Gz - dz*sum1 + tot*dzc*sum2
      end if
c-------------------------

c---------------------------------------
c     test0 = test0 + sum0   ! testing invariance with respect to ewald
c     test1 = test1 - 0.5*dzs*sum1 
c     test2 = test2 + toe*dzq*sum2
c---------------------------------------

c-------------------------------
      end if ! OF SUMMATION IN REAL SPACE
c-------------------------------

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

      p  = 0.0D0
      px = 0.0D0
      py = 0.0D0
      pz = 0.0D0

      Do i1=-Max2,Max2
       Do i2=-Max2,Max2

        k1 = i1*b11 + i2*b21 
        k2 = i1*b12 + i2*b22
 
        ks = k1*k1 + k2*k2
        k  = Dsqrt(ks)

        if(k.gt.0.000001) then  ! skip the zero point

c-------------------------
      if(method.eq.1) then          ! straight Fourier series
c-------------------------

        rho = Dabs(dzz)*k

        A = Dexp(-rho)/k

        if(Iopt.ne.1) then     ! for grad(G)
         B =  k1*A
         C =  k2*A
         D = -k *A
        end if

c-------------------------
      else if(method.eq.2) then     ! accelerated
c-------------------------

        w = k/(2.0D0*ew)
        T = 1.0D0/(1.0D0+0.5D0*w)          ! complementary error function
        ERFCC=    
     +  T*EXP(-W*W-1.26551223+T*(1.00002368+T*(.37409196+
     +  T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+
     +  T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
 
        fq0 =  erfcc/k 
        fq1 = -ks*fq0
        fq2 = -ks*fq1/9.0D0

        A = fq0 - 0.5D0 * dzs*fq1 + toe*dzq*fq2

        if(Iopt.ne.1) then    ! for grad(G)
         B = k1*A
         C = k2*A
         D =  - dz*fq1 + tot*dzc*fq2
        end if

c-----------------------------------
c       testt0 = testt0 + fq0
c       testt1 = testt1 - 0.5*dzs*fq1 
c       testt2 = testt2 + toe*dzq*fq2
c-----------------------------------

c-----------
      end if
c-----------

      end if

      arg = k1*dxx+k2*dyy
      cs  = Dcos(arg)

      p = p + A*cs

       if(Iopt.ne.1) then         ! for grad(G)
        sn = sin(arg)
        px = px - B*sn
        py = py - C*sn
        pz = pz + D*cs
       end if

c---------------------
c       p0 = p0+fq0*cs
c       p1 = p1+fq1*cs
c       p2 = p2+fq2*cs
c---------------------

c---
c end of computing grad(G)
c---

       End Do
      End Do

c-------------------------------
c END OF SUM IN RECIPROCAL SPACE
c-------------------------------

c----------
c FINISH UP
c----------

c--------------------
c pure Fourier series
c-------------------

      if(method.eq.1) then

        fc = 0.5/area
        G = fc*p
        G = G - fc*abs(dzz)     ! regularize by adding a linear function

        if(Iopt.ne.1) then      ! for grad(G)
         Gx = fc*px             ! normalize
         Gy = fc*py
         if(dzz.lt.0) pz = -pz
         Gz = fc*pz - fc*dzz/abs(dzz) ! regularize 
        end if

      end if

c-----------------
c expedited method
c-----------------

      if(method.eq.2) then

       prf = pi2/area

       p = p - 1.0D0/(ew*srpi)  ! contribution from zero wave number
       G = G + prf*p            ! add the sums in real and reciprocal space

       G = G/pi4              ! normalize

       if(Iopt.ne.1) then     ! for grad(G)

        Gx = Gx + prf*px      ! add the sums in real and reciprocal space
        Gy = Gy + prf*py
        Gz = Gz + prf*pz

        Gx = Gx/pi4           ! normalize
        Gy = Gy/pi4
        Gz = Gz/pi4

       end if

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

      end if

c-----
c Done
c-----

  100 Format (9(1x,f10.6))

      Return
      End
