      subroutine ft_nu 
     +
     +  (L
     +  ,k
     +  ,Npts
     +  ,x,y
     +  ,M
     +  ,af,bf
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

c------------------------------------------
c Computation of the complex Fourier coefficients
c of a periodic function.
c
c Integrals computed by the Simpson rule
c applied to overlapping parabolas
c
c SYMBOLS:
c -------
c
c L:    period
c k:    wave number = 2 pi /L
c npts: number of points
c x,y:	point coordinates
c
c M:    number of coefficients to be computed
c af,bf:   cosine and sine Fourier coefficients
c
c------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision L,k

      Integer p

      Dimension x(200),y(200)
      Dimension af(0:200),bf(0:200)

c-------------------------------
c Even/Odd number of intervals ?
c-------------------------------


      if(mod(Npts,2).eq.0) then

        write (6,*)
        write (6,*) "fourier_nu: ",npts," intervals: Even"

        M  = Npts/2
        cf = 0.5D0

      else

        write (6,*)
        write (6,*) "fourier_nu: ",npts," intervals: Odd"

        M  = (Npts-1)/2
        cf = 1.0D0

      end if

c-------
c extend
c-------

      x(Npts+1) = x(1)+L
      y(Npts+1) = y(1)
      x(Npts+2) = x(2)+L
      y(Npts+2) = y(2)

c----------
c launching
c----------

      Do p=0,M

       af(p) = 0.0D0
       bf(p) = 0.0D0

       Do i=1,Npts

c---
c complex transform
c by the Simpson rule
c---
        i1 = i+1
        i2 = i+2

        xhat = x(i)-x(1)
        arg  = p*k*xhat
        cs   = Dcos(arg)
        sn   = Dsin(arg)
        ff0  = cs*y(i)
        gg0  = sn*y(i)
        xx0  = xhat

        xhat = x(i1)-x(1)
        arg  = p*k*xhat
        cs   = Dcos(arg)
        sn   = Dsin(arg)
        ff1  = cs*y(i1)
        gg1  = sn*y(i1)
        xx1  = xhat

        xhat = x(i2)-x(1)
        arg  = p*k*xhat
        cs   = Dcos(arg)
        sn   = Dsin(arg)
        ff2  = cs*y(i2)
        gg2  = sn*y(i2)
        xx2  = xhat

        tm2 = (ff2-ff1)/(xx2-xx1)
        tm1 = (ff0-ff1)/(xx0-xx1)
        bbb = (tm2-tm1)/(xx2-xx0)
        ccc =  tm2-bbb*(xx2-xx1)

        dx2 = xx2-xx1
        dx0 = xx0-xx1

        af(p) = af(p) + ff1*(dx2   -dx0)
     +                + bbb*(dx2**3-dx0**3)/3.0D0
     +                + ccc*(dx2**2-dx0**2)/2.0D0

        tm2 = (gg2-gg1)/(xx2-xx1)
        tm1 = (gg0-gg1)/(xx0-xx1)
        bbb = (tm2-tm1)/(xx2-xx0)
        ccc =  tm2-bbb*(xx2-xx1)
        bf(p) = bf(p) + gg1*(dx2   -dx0)
     +                + bbb*(dx2**3-dx0**3)/3.0D0
     +                + ccc*(dx2**2-dx0**2)/2.0D0

       End Do

       af(p) = 0.5D0*af(p)       ! to account for overlap
       bf(p) = 0.5D0*bf(p)       ! to account for overlap

       af(p) = 2.0D0*af(p)/L
       bf(p) = 2.0D0*bf(p)/L

      End Do

c----------------------------
c Last coefficient is special
c----------------------------

      af(M) = cf*af(M)

c-----
c Done
c-----

      Return
      End
