      subroutine leg_fnc_ortho
     +
     +   (x
     +   ,j_max
     +   ,t
     +   )
      
c===========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

c--------------------------------------
c Evaluate the orthonormal set of 
c modified Legendre functions
c of degree j and order m
c
c for j=0,...,j_max, and m=0,...,j
c
c When multiplied by exp(i m phi),
c these modified Legendre functions
c comprise an orthonormal set over the surface
c of the unit sphere
c with unit weighting function,
c where i is the imaginary unit,
c phi is the azimuthal angle,
c x = cos(theta),
c theta is the meridional angle
c
c SYMBOLS
c -------
c
c x:       argument of the Legendre functions
c
c t(j,m):  Lower triangular matrix 
c          with dimensions (j_max) x (j_max)
c          holding the Legendre functions
c
c j_max:   maximum degree j
c
c--------------------------------------
      
      Implicit Double Precision(a-h,o-z)
      
      Dimension t(0:32,0:32)
      
      Parameter (tol=0.99999999D0)   ! for overflow

c----------
c constants
c----------

      pi = 3.14159 265358 D0

      half  = 0.5D0
      sqrth = Dsqrt(half)

c------------
c  initialize
c------------

      Do j=0,32
       Do m=0,32
         t(j,m) = 0.0D0
       End Do
      End Do
      
c---------
c  prepare
c---------

      if(x.gt. tol) x = tol  ! to prevent overflow
      if(x.lt.-tol) x =-tol  ! to prevent overflow

      xs = x*x
      sqrtxsc = dsqrt(1.0D0-xs)
      
c---------------------------
c Evaluate the matrix t
c column-by-column
c---------------------------

      Do m=0,j_max
        Do j=m,j_max+1
            
          r2jp1 = 2.0D0*j+1.0D0
          arg   = r2jp1/pi
          cf    = 0.5D0*sqrt(arg)

c---
c manual settings
c---

          if(m.eq.0.and.j.eq.0) then
             t(j,m) = cf
          else if(m.eq.0.and.j.eq.1) then
             t(j,m) = cf*x
          else if(m.eq.1.and.j.eq.1) then
             t(j,m) = cf*sqrth*sqrtxsc

c---
c recursion
c---
               
          else

            jpm = j+m
            jmm = j-m

            if(m.ne.j) then      ! use (A.7)

              arg = r2jp1*(2.0D0*j-1.0D0)*jmm/jpm
              c1  = sqrt(arg)
              arg = r2jp1*(jpm-1.0D0)*jmm*(jmm-1.0D0)
     +                                  /((2.0D0*j-3.0D0)*jpm)
              c2  = sqrt(arg)
              den = jmm-1.0D0+1.0D0

              t(j,m) = (c1*x*t(j-1,m)-c2*t(j-2,m))/den

            else      ! use (A.8)
                  
               arg = r2jp1*(jpm-1.0D0)*jpm/(2.0D0*j-1.0D0)
               c1  = sqrt(arg)
               arg = r2jp1*(jmm+1.0D0)*(jmm+2.0D0)/(2.0D0*j+3.0D0)
               c2  = sqrt(arg)
               den = r2jp1*sqrtxsc

               t(j,m) = (c1*t(j-1,m-1) - c2*t(j+1,m-1))/den
               
            end if                  

          end if

        End Do
      End Do
      
c---
c Display
c---
      
c     Do j=0,j_max
c       write (6,100) x,(t(j,m),m=0,j)
c     End Do
      
c-----
c Done
c-----
      
 100  Format (100(1x,f10.5))
   
      return
      end
