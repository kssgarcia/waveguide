      subroutine fft 
     +
     + (s
     + ,g
     + ,a,b
     + )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c           C. Pozrikidis
c Numerical Computation in Science and Engineering
c       Oxford University Press
c------------------------------------------------

c----------------------------------------
c Fast Fourier Transform on a (complex) 
c data set g(i).
c
c This subroutine returns the (complex) 
c frequency space representation.
c 
c Number of complex members of the data set g
c must be an integral power of 2:
c
c   n = 2^s
c
c Input data set is numbered from 0 to n-1
c
c The Fourier coefficients a and b are computed with
c the convention that elements 1 to n/2-1 
c correspond to positive frequency,
c and elements n/2+1 to n-1 correspond to negative
c frequency.  
c
c Element n/2 represents the nyquist frequency, both
c positive and negative.
c
c SYMBOLS:
c --------
c
c  g .... complex data set (g(i,1)=real, g(i,2)=imag.)
c  gr ... real part of frequency space representation
c  gi ... imaginary part of frequency space representation
c  s  ... size of data set (power of 2)
c  p .... ordering vector
c  mu ... integer power of two  
c  a .... real (cos) elements of frequency representation
c  b .... imaginary (sin) elements of frequency representation
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension g(0:4096,2),gr(0:4096),gsr(0:4096)
      Dimension             gi(0:4096),gsi(0:4096)
      Dimension a(0:4096),b(0:4096),cs(0:4096),sn(0:4096)

      Integer dft(0:4096,0:4096) ! optional

      Integer s
      Integer p(0:4096)

c----------
c constants
c----------

      pi = 3.14159 265358D0
      pi2 = 2.0D0*pi

c--------
c prepare
c--------

      n = 2**s    ! set size

c----------------------------
c copy the complex data set g 
c into the vectors gr and gi
c----------------------------

      Do i=0,n-1
       gr(i) = g(i,1)
       gi(i) = g(i,2)
      End Do

c--------------------------------
c Establish the ordering vector p 
c
c Algorithm (8.7.15)
c--------------------------------

      Do i=0,n-1
       p(i)=0
      End Do

      Do i=1,s
   
        mu = 2**(i-1)

        Do j=1,mu
         p(j) = 2*p(j)
        End Do

        Do j=0,mu-1
         p(j+mu) = p(j)+1
        End Do

c     write (6,*) (p(j),j=0,n-1)

      End Do

c------------------------------------------
c compute the sine and cosine of the angles
c------------------------------------------

      fc = pi2/n

      Do i=0,n-1
        tmp  = fc*i
        cs(i) = cos(tmp)
        sn(i) = sin(tmp) 
      End Do

c----------------
c perform the FFT 
c
c Algorithm 8.7.1
c----------------

      nsets = 1
      ndel  = n/2

      Do nstage=1,s

c---- ! optional
c      Do i=0,n-1 
c      Do j=0,n-1
c       dft(i,j)=111
c      End Do
c      End Do
c---
   
       ig = 0
       Do iset=1,nsets
        nrt = n/nsets-1
        Do m=0,nrt
         j = Mod(m,ndel) + 2*ndel*(iset-1)
         l = int(ig/ndel)
         l = p(l)
c--- optional
c        dft(ig,j)=0
c        dft(ig,j+ndel)=l
c        write (6,*) j+1,j+ndel+1,l
c---
         gsr(ig) = gr(j)+cs(l)*gr(j+ndel)-sn(l)*gi(j+ndel)
         gsi(ig) = gi(j)+cs(l)*gi(j+ndel)+sn(l)*gr(j+ndel)
         ig = ig+1
         End Do
        End Do

c--- optional
c       write (6,*)
c       Do i=0,n-1 
c        write (6,*) (dft(i,j),j=0,n-1)
c       End Do
c--- 

        Do i=0,n-1
          gr(i) = gsr(i)
          gi(i) = gsi(i)
        End Do
 
        nsets = 2*nsets
        ndel  = ndel/2

      End Do

c     stop

c--------------------------------------
c unscramble with the indexing vector p
c--------------------------------------

      Do i=0,n-1
       a(i) = 2.0D0*gr(p(i))
       b(i) = 2.0D0*gi(p(i))
      End Do

c-------
c adjust
c-------

      a(n/2)=0.5D0*a(n/2);
      b(n/2)=0.5D0*b(n/2);

c-----
c Done
c-----

      Return
      End
