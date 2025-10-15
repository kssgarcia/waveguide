      program binary

c===========================================
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

c------------------------------------------------
c This program accompanies the book:
c           C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c        Oxford University Press
c------------------------------------------------

c-----------------------------------
c Convert a positive real number to the
c corresponding binary number
c
c SYMBOLS:
c --------
c
c nb:   array of binary digits
c v:    auxiliary variable
c k:    number of binary digits
c       corrresponding to the integral part
c l:    number of binary digits
c       corrresponding to the decimal part
c---------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension v(128),nb(-16:16)

c----------
c constants
c----------

      two = 2.0D0

c-------
c launch
c-------

  96  Continue

      write (6,*)
      write (6,*) " Will convert a positive real number"
      write (6,*) " into a binary"
      write (6,*) " Enter the number (e.g., 10.45)"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read (5,*) a

      if(a.eq.0) Go to 99    ! quit

c---
c initialize
c---

      Do i=-16,16
       nb(i) = 0;
      End Do

c---------------------
c extract the integral
c and decimal parts
c---------------------

      aint = int(a)
      adec = a-aint

c-------------------------
c handle the integral part
c-------------------------

      Do i=1,32
       v(i)    = mod(aint,two)
       nb(i-1) = nint(v(i))
       aint    = 0.5D0*(aint-v(i))
       if(aint.le.0.00001) Go to 97
      End Do

 97   Continue

      k = i-1    ! number of binary digits corresponding
                 ! to the integral part
c------------------------
c handle the decimal part
c------------------------

      adec = 2.0*adec

      Do i=1,50
       v(i)   = int(adec)
       nb(-i) = v(i)
       adec   = 2.0D0*(adec-v(i))
       if(adec.le.0.0000001) go to 98
      End Do

 98   Continue

      l = i      ! number of binary digits corresponding
                 ! to the decimal part
c------------------------
c print the binary number
c------------------------

      write (6,*) 
      write (6,*) " The binary form is:"
      write (6,*) 
      write (6,100) (nb(j),j= 15, 0,-1),(nb(j),j=-1,-16,-1)

c-----------------------------
c reproduce the decimal number
c by Horner's algorithm
c-----------------------------

      c = nb(k)

      Do i=k-1,0,-1
       c = 2.0D0*c+nb(i)
      End Do

      d = nb(-l)

      Do i=-l+1,-1
       d = 0.5D0*d+nb(i)
      End Do

      e = c+0.5D0*d

      write (6,*) 
      write (6,*) "Original and reproduced numbers:"
      write (6,*) 
      write (6,*) a,e

c-----------------
c return to repeat
c-----------------

      Go to 96   ! repeat

c-----
c done
c-----

 99   Continue

100   Format (16(i1),1x,16(i1))

      stop
      end
