      program Fourier_ident

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
c------------------------------------------------

c---------------------------------------------
c Verify the fourier identities of Table 8.6.2
c---------------------------------------------

      Implicit Double Precision (a-h, o-z)

      Double Precision L,k
      Integer p,q

c----------
c constants
c----------

      pi  = 3.14159 265358
      pi2 = 2.0*pi

c---
c input
c---

      write (6,*)
      write (6,*) " Please enter the end-points a, b"
      write (6,*) " --------------------------------"
      read  (5,*) a, b
      write (6,*)

      write (6,*) 
      write (6,*) " Please enter the number of intervals N"
      write (6,*) "                              0 to quit"
      write (6,*) " --------------------------------------"
      read  (5,*) N

      If(N.eq.0) Go to 99

 98   Continue

      write (6,*)
      write (6,*) " Please enter p and q (both integers)"
      write (6,*) "                           99 to quit"
      write (6,*) " ------------------------------------"
      read (5,*) p,q

      If(p.eq.99) Go to 99
      If(p.eq.99) Go to 99

c---
c prepare
c---

      L = b-a
      k = pi2/L
      h = L/(N-1.0+1.0)

c---
c five identities
c---

      sumcp = 0.0
      sumsp = 0.0
      sumcq = 0.0
      sumsq = 0.0
      sumcc = 0.0
      sumss = 0.0
      sumcs = 0.0

      Do j = 1,N

        xhj  = (j-1)*h
        argp = p*k*xhj
        argq = q*k*xhj

        cp    = cos(argp)
        sp    = sin(argp)
        cq    = cos(argq)
        sq    = sin(argq)
        sumcp = sumcp + cp
        sumsp = sumsp + sp
        sumcq = sumcq + cq
        sumsq = sumsq + sq
        sumcc = sumcc + cp*cq
        sumss = sumss + sp*sq
        sumcs = sumcs + cp*sq

      End Do

c---
c printing
c---

      write (6,110) sumcp
      write (6,111) sumsp
      write (6,112) sumcq
      write (6,113) sumsq
      write (6,114) sumcc
      write (6,115) sumss
      write (6,116) sumcs

      Go to 98   ! repeat

c---
c Done
c---

 99   Continue

 100  Format (1x,f10.5,1x,f10.5,1x,f10.5)
 110  Format ( "Sum of cosines for p ", f10.5)
 111  Format ( "Sum of   sines for p ", f10.5)
 112  Format ( "Sum of cosines for q ", f10.5)
 113  Format ( "Sum of   sines for q ", f10.5)
 114  Format ( "Sum of cos-cos       ", f10.5)
 115  Format ( "Sum of sin-sin       ", f10.5)
 116  Format ( "Sum of cos-sin       ", f10.5)

      Stop
      End
