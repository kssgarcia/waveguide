      program spec_fnc

c===========================================
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

c------------------------------------------------
c This program accompanies the book:
c            C. Pozrikidis
c Numerical Computation in Science and Engineering
c         Oxford University Press
c------------------------------------------------

c-----------------------------------------
c  Draw profiles of special functions
c-----------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      null = 0

c--------
c prepare
c--------

      open (1,file="spec_fnc.out")

c------------
c preferences
c------------

  98  Continue

      write (6,*)
      write (6,*) "  Enter:"
      write (6,*)
      write (6,*) "  1 for the complete elliptic integral"
      write (6,*) "    of the first kind"
      write (6,*) "  2 for the complete elliptic integral"
      write (6,*) "    of the second kind"
      write (6,*) "  3 for the error function"
      write (6,*) "  4 for the complementary error function"
      write (6,*) "  5 for the Bessel function J0"
      write (6,*) "  6 for the Bessel function J1"
      write (6,*) "  7 for the Kelvin function ber_0"
      write (6,*) "  8 for the Kelvin function bei_0"
      write (6,*) "  9 for the modified Bessel function I_0"
      write (6,*) " 10 for the modified Bessel function I_1"
      write (6,*) " 11 for the modified Bessel function K_0"
      write (6,*) " 12 for the modified Bessel function K_1"
      write (6,*) " 0  to quit"
      write (6,*) " ----------"
      read  (5,*) menu

      if(menu.eq.0) Go to 99

c-----------------------
      if(menu.eq.1) then        ! complete elliptic integral F
c-----------------------

       Nprof  = 128
       Nprof1 = Nprof+1
       xmin   = 0.001
       xmax   = 0.999
       step   = (xmax-xmin)/Nprof

       write (6,100) Nprof1
       write (1,100) Nprof1

       Do i=1,Nprof1
        x  = xmin+(i-1.0)*step
        x2 = x**2
        call ell_int (x2,F,E)
        write (6,100) i,x2,F
        write (1,100) i,x2,F
       End Do

c----------------------------
      else if(menu.eq.2) then       ! complete elliptic integral F
c----------------------------

       Nprof  = 128
       Nprof1 = Nprof+1
       xmin   = 0.001
       xmax   = 0.999
       step   = (xmax-xmin)/Nprof

       write (6,100) Nprof1
       write (1,100) Nprof1

       Do i=1,Nprof1
        x  = xmin+(i-1.0)*step
        x2 = x**2
        call ell_int (x2,F,E)
        write (6,100) i,x2,E
        write (1,100) i,x2,E
       End Do

c----------------------------
      else if(menu.eq.3) then            ! error function
c----------------------------

       Nprof  = 128
       Nprof1 = Nprof+1
       xmin   = 0.001
       xmax   = 2.0
       step   = (xmax-xmin)/Nprof

       write (6,100) Nprof1
       write (1,100) Nprof1

       Do i=1,Nprof1
        x = xmin+(i-1.0D0)*step
        y = erfun(x)
        write (6,100) i,x,y
        write (1,100) i,x,y
       End Do

c----------------------------
      else if(menu.eq.4) then  ! complementary error function
c----------------------------

       Nprof  = 128
       Nprof1 = Nprof+1
       xmin   = 0.001
       xmax   = 2.0
       step   = (xmax-xmin)/Nprof

       write (6,100) Nprof1
       write (1,100) Nprof1

       Do i=1,Nprof1
        x = xmin+(i-1.0D0)*step
        y = 1.0D0-erfun(x)
        write (6,100) i,x,y
        write (1,100) i,x,y
       End Do

c----------------------------
      else if(menu.eq.5) then         ! Bessel function J0
c----------------------------

       Nprof  = 128
       Nprof1 = Nprof+1
       xmin   = 0.001
       xmax   = 16.0
       step   = (xmax-xmin)/Nprof

       write (6,100) Nprof1
       write (1,100) Nprof1

       Do i=1,Nprof1
        x = xmin+(i-1.0)*step
        y = bess_J0(x)
        write (6,100) i,x,y
        write (1,100) i,x,y
       End Do

c----------------------------
      else if(menu.eq.6) then               ! Bessel function J1
c----------------------------

       Nprof  = 128
       Nprof1 = Nprof+1
       xmin   = 0.001
       xmax   = 16.0
       step   = (xmax-xmin)/nprof

       write (1,100) Nprof1
       write (6,100) Nprof1

       Do i=1,Nprof1
        x = xmin+(i-1.0)*step
        y = bess_J1(x)
        write (6,100) i,x,y
        write (1,100) i,x,y
       End Do

c----------------------------
      else if(menu.eq.7) then  ! Kelvin function ber_0
c----------------------------

       Nprof  = 64
       Nprof1 = Nprof+1
       xmin   = 0.001
       xmax   = 8.0
       step   = (xmax-xmin)/Nprof

       write (6,100) Nprof1
       write (1,100) Nprof1

       Do i=1,Nprof1

        x = xmin+(i-1.0)*step

        call ber_bei_0
     +
     +    (Iopt
     +    ,X
     +    ,ber_0,bei_0
     +    ,ber_0_p,bei_0_p
     +    )

        y = ber_0

        write (6,100) i,x,y
        write (1,100) i,x,y

       End Do

c----------------------------
      else if(menu.eq.8) then         ! Kelvin function bei_0
c----------------------------

       Nprof  = 64
       Nprof1 = Nprof+1
       xmin   = 0.001
       xmax   = 8.0
       step   = (xmax-xmin)/nprof

       write (6,100) Nprof1
       write (1,100) Nprof1

       Do i=1,Nprof1

        x = xmin+(i-1.0)*step

        call ber_bei_0
     +
     +    (Iopt
     +    ,X
     +    ,ber_0,bei_0
     +    ,ber_0_p,bei_0_p
     +    )

        y = bei_0

        write (6,100) i,x,y
        write (1,100) i,x,y

       End Do

c----------------------------
      else if(menu.eq.9) then    ! modified Bessel function I_0
c----------------------------

       Nprof  = 64
       Nprof1 = Nprof+1
       xmin   = 0.001
       xmax   = 8.0
       step   = (xmax-xmin)/nprof

       write (6,100) Nprof1
       write (1,100) Nprof1

       Do i=1,Nprof1

        x = xmin+(i-1.0)*step

        call bess_I01K01
     +
     +     (x
     +     ,Iswitch
     +     ,BI0,BI1
     +     ,BK0,BK1
     +     )

        y = BI0

        write (6,100) i,x,y
        write (1,100) i,x,y

       End Do

c-----------------------------
      else if(menu.eq.10) then      ! modified Bessel function I_1
c-----------------------------

       Nprof  = 64
       Nprof1 = Nprof+1
       xmin   = 0.001
       xmax   = 8.0
       step   = (xmax-xmin)/nprof

       write (6,100) Nprof1
       write (1,100) Nprof1

       Do i=1,Nprof1

        x = xmin+(i-1.0)*step

        call bess_I01K01
     +
     +     (x
     +     ,Iswitch
     +     ,BI0,BI1
     +     ,BK0,BK1
     +     )

        y = BI1

        write (6,100) i,x,y
        write (1,100) i,x,y
       End Do

c----------------------------
      else if(menu.eq.11) then         ! modified Bessel function K_0
c----------------------------

       Nprof  = 64
       Nprof1 = Nprof+1
       xmin   = 0.001
       xmax   = 8.0
       step   = (xmax-xmin)/nprof

       write (6,100) Nprof1
       write (1,100) Nprof1

       Do i=1,Nprof1
        x = xmin+(i-1.0)*step
        call bess_I01K01
     +
     +    (x
     +    ,Iswitch
     +    ,BI0,BI1
     +    ,BK0,BK1
     +    )
        y = BK0
        write (6,100) i,x,y
        write (1,100) i,x,y
       End Do

c-----------------------------
      else if(menu.eq.12) then      ! modified Bessel function K_1
c-----------------------------

       Nprof  = 64
       Nprof1 = Nprof+1
       xmin   = 0.001
       xmax   = 8.0
       step   = (xmax-xmin)/nprof

       write (6,100) Nprof1
       write (1,100) Nprof1

       Do i=1,Nprof1

        x = xmin+(i-1.0)*step

        call bess_I01K01
     +
     +    (x
     +    ,Iswitch
     +    ,BI0,BI1
     +    ,BK0,BK1
     +    )

        y = BK1

        write (6,100) i,x,y
        write (1,100) i,x,y

       End Do

c-----------
      End If
c-----------

      Go to 98   ! repeat

c-----
c done
c-----

  99  Continue

      write (1,*) null
      close (1)

 100  Format (i5,10(1x,f15.10))

      stop
      end
