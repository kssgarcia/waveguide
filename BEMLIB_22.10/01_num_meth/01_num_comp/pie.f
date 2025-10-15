      program pie

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
c          C. Pozrikidis
c Numerical Computation in Science and Engineering
c         Oxford University Press
c------------------------------------------------

c--------------------------------------------------
c Compute the number pi as the limit 
c of a sequence of the perimeter of a regular
c n-sided ploygons inscribed within the unit circle,
c in single and double precision
c
c The recurrence relation between succesive elemnts 
c of the series is given by:
c
c              m+1/2            2    2m+2
c       p   = 2      (1-sqrt(1-p  / 2     ))
c        m+1                    m
c------------------------------------------------

      Double Precision argd,twod,srtd,pd,pid

      Character*20 filename

c------------
c output file
c------------

c     write (6,10)
c     read  (5,'(A20)') filename
c     open (1,file=filename,status='unknown')

c----------
c constants
c----------

      two  = 2.0               ! single precision
      srt  = sqrt(two)

      twod = 2.0D0             ! double precision
      srtd = Dsqrt(twod)

c------------
c exact value
c------------

      pi  = 3.141592            ! single precision
      pid = 3.141592 65358D0    ! double precision

      write (6,*)
      write (6,*) " Exact values in single and double precision:"
      write (6,*)
      write (6,101) pi,pid

 98   Continue

      write (6,*)
      write (6,*) "How many iterations? (0 to quit)"
      write (6,*) "--------------------------------"
      read  (5,*) N

      If(N.eq.0) Go to 99

c-----------------
c single precision
c-----------------

      write (6,*)
      write (6,*) " Sequence in single precision"
      write (6,*)

      p = 2.0*srt
      m = 1
      write (6,100) m,p
c     write (1,100) m,p

      Do m=2,N
          arg = 1.0-(p/2.0**m)**2
        If(arg.lt.0.0) Go to 91
          tmp = sqrt(arg)
        If(tmp.gt.1.0) Go to 91
          p = srt  * 2.0**m * sqrt(1.0-tmp)
        write (6,100) m,p
c       write (1,100) m,p
      End Do

  91  Continue

c-----------------
c double precision
c-----------------

      write (6,*)
      write (6,*) " Sequence in double precision:"
      write (6,*)

      pd = 2.0D0*srtd

      m = 1
      write (6,100) m,pd
c     write (1,100) m,pd

      Do m=2,N
        argd = 1.0D0-(pd/2.0D0**m)**2
        If(argd.lt.0.0) Go to 92
        tmpd = Dsqrt(argd)
        If(tmpd.gt.1.0D0) Go to 92
        pd  = srtd * 2.0D0**m * dsqrt(1.0D0-tmpd)
        write (6,100) m,pd
c       write (1,100) m,pd
      End Do

 92   Continue

c---
c Return to repeat
c---

      Go to 98

c-----
c Done
c-----

  99  Continue

c     close(1)

 10   Format (' Please enter the name of output file: ',$)
 100  Format (1x,i4,1x,f16.14,1x,f16.14)
 101  Format (1x,f16.14,5x,f16.14)

      Stop
      End
