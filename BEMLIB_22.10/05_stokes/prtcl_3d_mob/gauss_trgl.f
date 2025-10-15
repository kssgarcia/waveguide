      subroutine trgl_quad (n,xi,eta,w)

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------
c Abscissas and Weight Factors for 
c Gaussian Integration over a triangle
c------------------------------------

      Implicit double precision (a-h,o-z)
      Double precision kappa,mu

      Dimension xi(20),eta(20),w(20)

c------------
c Preferences
c------------

      If(   n.ne.  1
     +  .and.n.ne. 3
     +  .and.n.ne. 4
     +  .and.n.ne. 6
     +  .and.n.ne. 7
     +  .and.n.ne. 9
     +  .and.n.ne.12
     +  .and.n.ne.13) then
         write (6,*)
         write (6,*) ' Selected number of triangular gaussian',
     +               ' points not available'
         write (6,*) ' will take n=7'
         write (6,*)
         n=7
      End If      

c-----------
      If(n.eq.1) then
c-----------

         xi (1) = 1.0/3.0
         eta(1) = 1.0/3.0
         w  (1) = 1.0

c-----------
      Else If(n.eq.3) then
c-----------

         xi (1) = 1.0/6.0
         xi (2) = 4.0/6.0
         xi (3) = 1.0/6.0
         eta(1) = 1.0/6.0
         eta(2) = 1.0/6.0
         eta(3) = 4.0/6.0
         w  (1) = 1.0/3.0
         w  (2) = 1.0/3.0
         w  (3) = 1.0/3.0

c-----------
      Else If(n.eq.4) then
c-----------

         xi(1)  =  1.0/ 5.0
         xi(2)  =  3.0/ 5.0
         xi(3)  =  1.0/ 5.0
         xi(4)  =  1.0/ 3.0
         eta(1) =  1.0/ 5.0
         eta(2) =  1.0/ 5.0
         eta(3) =  3.0/ 5.0
         eta(4) =  1.0/ 3.0
         w(1)   = 25.0/48.0
         w(2)   = 25.0/48.0
         w(3)   = 25.0/48.0
         w(4)   =-27.0/48.0

c-----------
      Else If(n.eq.6) then
c-----------

         alpha  = 0.81684 75729 80459
         beta   = 0.44594 84909 15965
         gamma  = 0.10810 30181 68070
         delta  = 0.09157 62135 09771
         omega1 = 0.10995 17436 55322
         omega2 = 0.22338 15896 78011

         xi(1)  = delta
         xi(2)  = alpha
         xi(3)  = delta
         xi(4)  = beta
         xi(5)  = gamma
         xi(6)  = beta
         eta(1) = delta
         eta(2) = delta
         eta(3) = alpha
         eta(4) = beta
         eta(5) = beta
         eta(6) = gamma
         w(1)   = omega1
         w(2)   = omega1
         w(3)   = omega1
         w(4)   = omega2
         w(5)   = omega2
         w(6)   = omega2

c-----------
      Else If(n.eq.7) then
c-----------

         alpha  = 0.79742 69583 53087
         beta   = 0.47014 20641 05115
         gamma  = 0.05971 58717 89770
         delta  = 0.10128 65073 23456
         omega1 = 0.12593 91805 44827
         omega2 = 0.13239 41527 88506

         xi(1)  = delta 
         xi(2)  = alpha
         xi(3)  = delta
         xi(4)  = beta
         xi(5)  = gamma
         xi(6)  = beta
         xi(7)  = 1.0/3.0
         eta(1) = delta
         eta(2) = delta
         eta(3) = alpha
         eta(4) = beta
         eta(5) = beta
         eta(6) = gamma
         eta(7) = 1.0/3.0
         w(1)   = omega1
         w(2)   = omega1
         w(3)   = omega1
         w(4)   = omega2
         w(5)   = omega2
         w(6)   = omega2
         w(7)   = 0.2250

c-----------
      Else If(n.eq.9) then
c-----------

         alpha  = 0.12494 95032 33232
         kappa  = 0.16540 99273 89841
         rho    = 0.79711 26518 60071
         delta  = 0.43752 52483 83384
         mu     = 0.03747 74207 50088
         omega1 = 0.20595 05047 60887
         omega2 = 0.06369 14142 86223
         xi(1)  = delta
         xi(2)  = alpha
         xi(3)  = delta
         xi(4)  = kappa
         xi(5)  = mu
         xi(6)  = rho
         xi(7)  = kappa
         xi(8)  = mu
         xi(9)  = rho
         eta(1) = delta
         eta(2) = delta
         eta(3) = alpha
         eta(4) = mu
         eta(5) = kappa
         eta(6) = kappa
         eta(7) = rho
         eta(8) = rho
         eta(9) = mu
         w(1)   = omega1
         w(2)   = omega1
         w(3)   = omega1
         w(4)   = omega2
         w(5)   = omega2
         w(6)   = omega2
         w(7)   = omega2
         w(8)   = omega2
         w(9)   = omega2

c-----------
      Else If(n.eq.12) then
c-----------

         alpha  = 0.87382 19710 16996
         beta   = 0.24928 67451 70910
         gamma  = 0.50142 65096 58179
         delta  = 0.06308 90144 91502
         rho    = 0.63650 24991 21399
         kappa  = 0.31035 24510 33785
         mu     = 0.05314 50498 44816
         omega1 = 0.05084 49063 70207
         omega2 = 0.11678 62757 26379
         omega3 = 0.08285 10756 18374

         xi(1)  = delta
         xi(2)  = alpha
         xi(3)  = delta
         xi(4)  = beta
         xi(5)  = gamma
         xi(6)  = beta
         xi(7)  = kappa
         xi(8)  = mu
         xi(9)  = rho
         xi(10) = kappa
         xi(11) = mu
         xi(12) = rho
         eta(1) = delta
         eta(2) = delta
         eta(3) = alpha
         eta(4) = beta
         eta(5) = beta
         eta(6) = gamma
         eta(7) = mu
         eta(8) = kappa
         eta(9) = kappa
         eta(10)= rho
         eta(11)= rho
         eta(12)= mu
         w(1)   = omega1
         w(2)   = omega1
         w(3)   = omega1
         w(4)   = omega2
         w(5)   = omega2
         w(6)   = omega2
         w(7)   = omega3
         w(8)   = omega3
         w(9)   = omega3
         w(10)  = omega3
         w(11)  = omega3
         w(12)  = omega3

c-----------
      Else If (n.eq.13) then
c-----------

         alpha  = 0.47930 80678 41923
         beta   = 0.06513 01029 02216
         gamma  = 0.86973 97941 95568
         delta  = 0.26034 59660 79038
         rho    = 0.63844 41885 69809
         kappa  = 0.31286 54960 04875
         mu     = 0.04869 03154 25316
         omega1 = 0.17561 52574 33204
         omega2 = 0.05334 72356 08839
         omega3 = 0.07711 37608 90257
         omega4 =-0.14957 00444 67670
         xi(1)  = delta
         xi(2)  = alpha
         xi(3)  = delta
         xi(4)  = beta
         xi(5)  = gamma
         xi(6)  = beta
         xi(7)  = kappa
         xi(8)  = mu
         xi(9)  = rho
         xi(10) = kappa
         xi(11) = mu
         xi(12) = rho
         xi(13) = 1.d0/3.d0
         eta(1) = delta
         eta(2) = delta
         eta(3) = alpha
         eta(4) = beta
         eta(5) = beta
         eta(6) = gamma
         eta(7) = mu
         eta(8) = kappa
         eta(9) = kappa
         eta(10)= rho
         eta(11)= rho
         eta(12)= mu
         eta(13)= 1.d0/3.d0
         w(1)   = omega1
         w(2)   = omega1
         w(3)   = omega1
         w(4)   = omega2
         w(5)   = omega2
         w(6)   = omega2
         w(7)   = omega3
         w(8)   = omega3
         w(9)   = omega3
         w(10)  = omega3
         w(11)  = omega3
         w(12)  = omega3
         w(13)  = omega4

c-----------
      End If
c-----------

      Return
      End
