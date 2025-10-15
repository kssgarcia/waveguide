      subroutine gauss_trgl 
     +
     +  (N
     +  ,xi,eta,w
     +  )

c===========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

c------------------------------------
c Abscissas (xi, eta) and weights (w)
c for Gaussian integration over a 
c flat triangle
c
c The integration is performed with respect
c to the triangle barycentric coordinates
c
c SYMBOLS:
c -------
c
c N: order of the quadrature
c    choose from 1,3,4,5,7,9,12,13
c
c    Default value is 7
c------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension xi(20),eta(20),w(20)

c------
c traps
c------

      if(    N.ne. 1
     +  .and.N.ne. 3
     +  .and.N.ne. 4
     +  .and.N.ne. 6
     +  .and.N.ne. 7
     +  .and.N.ne. 9
     +  .and.N.ne.12
     +  .and.N.ne.13) then

        write (6,*)
        write (6,*) ' gauss_trgl:'
        write (6,*)
        write (6,*) ' Number of Guass triangle quadrature',
     +              ' points is not available; will take N=7'
        write (6,*)
        N = 7
      end if      

c--------------------
      if(N.eq.1) then
c--------------------

       xi(1) = 1.0D0/3.0D0
      eta(1) = 1.0D0/3.0D0
        w(1) = 1.0D0

c-------------------------
      else if(N.eq.3) then
c-------------------------

       xi(1) = 1.0D0/6.0D0
      eta(1) = 1.0D0/6.0D0
        w(1) = 1.0D0/3.0D0

       xi(2) = 2.0D0/3.0D0
      eta(2) = 1.0D0/6.0D0
        w(2) = w(1)

       xi(3) = 1.0D0/6.0D0
      eta(3) = 2.0D0/3.0D0
        w(3) = w(1)

c-------------------------
      else if(N.eq.4) then
c-------------------------

       xi(1) =   1.0D0/3.0D0
      eta(1) =   1.0D0/3.0D0
        w(1) = -27.0D0/48.0D0

       xi(2) =  1.0D0/5.0D0
      eta(2) =  1.0D0/5.0D0
        w(2) = 25.0D0/48.0D0

       xi(3) =  3.0D0/5.0D0
      eta(3) =  1.0D0/5.0D0
        w(3) = 25.0D0/48.0D0

       xi(4) =  1.0D0/5.0D0
      eta(4) =  3.0D0/5.0D0
        w(4) = 25.0D0/48.0D0

c-------------------------
      else if(N.eq.6) then
c-------------------------

      al = 0.81684 75729 80459 D0
      be = 0.44594 84909 15965 D0
      ga = 0.10810 30181 68070 D0
      de = 0.09157 62135 09771 D0
      o1 = 0.10995 17436 55322 D0
      o2 = 0.22338 15896 78011 D0

      xi(1) = de
      xi(2) = al
      xi(3) = de
      xi(4) = be
      xi(5) = ga
      xi(6) = be

      eta(1) = de
      eta(2) = de
      eta(3) = al
      eta(4) = be
      eta(5) = be
      eta(6) = ga

      w(1) = o1
      w(2) = o1
      w(3) = o1
      w(4) = o2
      w(5) = o2
      w(6) = o2

c-------------------------
      else if(N.eq.7) then
c-------------------------

      al = 0.79742 69583 53087 D0
      be = 0.47014 20641 05115 D0     
      ga = 0.05971 58717 89770 D0
      de = 0.10128 65073 23456 D0
      o1 = 0.12593 91805 44827 D0
      o2 = 0.13239 41527 88506 D0

      xi(1) = de
      xi(2) = al
      xi(3) = de
      xi(4) = be
      xi(5) = ga
      xi(6) = be
      xi(7) = 1.0D0/3.0D0

      eta(1) = de
      eta(2) = de
      eta(3) = al
      eta(4) = be
      eta(5) = be
      eta(6) = ga
      eta(7) = 1.0D0/3.0D0

      w(1) = o1
      w(2) = o1
      w(3) = o1
      w(4) = o2
      w(5) = o2
      w(6) = o2
      w(7) = 0.225 D0

c-------------------------
      else if(N.eq.9) then
c-------------------------

      al = 0.12494 95032 33232 D0
      qa = 0.16540 99273 89841 D0
      rh = 0.79711 26518 60071 D0
      de = 0.43752 52483 83384 D0
      ru = 0.03747 74207 50088 D0
      o1 = 0.20595 05047 60887 D0
      o2 = 0.06369 14142 86223 D0

      xi(1) = de
      xi(2) = al
      xi(3) = de
      xi(4) = qa
      xi(5) = ru
      xi(6) = rh
      xi(7) = qa
      xi(8) = ru
      xi(9) = rh

      eta(1) = de
      eta(2) = de
      eta(3) = al
      eta(4) = ru
      eta(5) = qa
      eta(6) = qa
      eta(7) = rh
      eta(8) = rh
      eta(9) = ru

      w(1) = o1
      w(2) = o1
      w(3) = o1
      w(4) = o2
      w(5) = o2
      w(6) = o2
      w(7) = o2
      w(8) = o2
      w(9) = o2

c--------------------------
      else if(N.eq.12) then
c--------------------------

      al = 0.87382 19710 16996 D0
      be = 0.24928 67451 70910 D0
      ga = 0.50142 65096 58179 D0
      de = 0.06308 90144 91502 D0
      rh = 0.63650 24991 21399 D0
      qa = 0.31035 24510 33785 D0
      ru = 0.05314 50498 44816 D0
      o1 = 0.05084 49063 70207 D0
      o2 = 0.11678 62757 26379 D0
      o3 = 0.08285 10756 18374 D0

      xi(1)  = de
      xi(2)  = al
      xi(3)  = de
      xi(4)  = be
      xi(5)  = ga
      xi(6)  = be
      xi(7)  = qa
      xi(8)  = ru
      xi(9)  = rh
      xi(10) = qa
      xi(11) = ru
      xi(12) = rh

      eta(1)  = de
      eta(2)  = de
      eta(3)  = al
      eta(4)  = be
      eta(5)  = be
      eta(6)  = ga
      eta(7)  = ru
      eta(8)  = qa
      eta(9)  = qa
      eta(10) = rh
      eta(11) = rh
      eta(12) = ru

      w(1)  = o1
      w(2)  = o1
      w(3)  = o1
      w(4)  = o2
      w(5)  = o2
      w(6)  = o2
      w(7)  = o3
      w(8)  = o3
      w(9)  = o3
      w(10) = o3
      w(11) = o3
      w(12) = o3

c---------------------------
      else if (N.eq.13) then
c---------------------------

      al = 0.47930 80678 41923 D0
      be = 0.06513 01029 02216 D0
      ga = 0.86973 97941 95568 D0
      de = 0.26034 59660 79038 D0
      rh = 0.63844 41885 69809 D0
      qa = 0.31286 54960 04875 D0
      ru = 0.04869 03154 25316 D0
      o1 = 0.17561 52574 33204 D0
      o2 = 0.05334 72356 08839 D0
      o3 = 0.07711 37608 90257 D0
      o4 =-0.14957 00444 67670 D0

      xi(1)  = de
      xi(2)  = al
      xi(3)  = de
      xi(4)  = be
      xi(5)  = ga
      xi(6)  = be
      xi(7)  = qa
      xi(8)  = ru
      xi(9)  = rh
      xi(10) = qa
      xi(11) = ru
      xi(12) = rh
      xi(13) = 1.0D0/3.0D0

      eta(1)  = de
      eta(2)  = de
      eta(3)  = al
      eta(4)  = be
      eta(5)  = be
      eta(6)  = ga
      eta(7)  = ru
      eta(8)  = qa
      eta(9)  = qa
      eta(10) = rh
      eta(11) = rh
      eta(12) = ru
      eta(13) = 1.0D0/3.0D0

      w(1)  = o1
      w(2)  = o1
      w(3)  = o1
      w(4)  = o2
      w(5)  = o2
      w(6)  = o2
      w(7)  = o3
      w(8)  = o3
      w(9)  = o3
      w(10) = o3
      w(11) = o3
      w(12) = o3
      w(13) = o4

c-----------
      end if
c-----------

      return
      end
