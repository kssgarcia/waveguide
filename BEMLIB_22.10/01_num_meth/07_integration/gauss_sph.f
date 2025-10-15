      subroutine gauss_sph 
     +
     +  (Nint
     +  ,xx,yy,zz
     +  ,ww
     +  )

c==============================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==============================================

c------------------------------------------------
c This program accompanies the book:
c
c                 C. Pozrikidis
c Numerical Computation in Science and Engineering
c           Oxford University Press
c------------------------------------------------

c----------------------------------------
c Abscissae and weights for integration 
c over the surface of a sphere 
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension xx(50),yy(50),zz(50),ww(50)

      if(Nint.ne.6.and.Nint.ne.18) then
         write (6,*)
         write (6,*) " gauss_sph: Requested number of points"
         write (6,*) "            is not available "
         write (6,*) "            Will use 18 points"
         Nint = 18
      end if

c-----------------------
      if(Nint.eq.6) then
c-----------------------

       xx(1) = 1.0D0
       yy(1) = 0.0D0
       zz(1) = 0.0D0
       ww(1) = 1.0D0/6.0D0

       xx(2) =-1.0D0
       yy(2) = 0.0D0
       zz(2) = 0.0D0
       ww(2) = ww(1)

       xx(3) = 0.0D0
       yy(3) = 1.0D0
       zz(3) = 0.0D0
       ww(3) = ww(1)

       xx(4) = 0.0D0
       yy(4) =-1.0D0
       zz(4) = 0.0D0
       ww(4) = ww(1)

       xx(5) = 0.0D0
       yy(5) = 0.0D0
       zz(5) = 1.0D0
       ww(5) = ww(1)

       xx(6) = 0.0D0
       yy(6) = 0.0D0
       zz(6) =-1.0D0
       ww(6) = ww(1)

c------------------------------
      else if (Nint.eq.18) then
c------------------------------

       two  = 2.0D0
       srti = 1.0D0/sqrt(two)

       xx(1) = 1.0D0
       yy(1) = 0.0D0
       zz(1) = 0.0D0
       ww(1) = 1.0D0/30.0D0

       xx(2) =-1.0D0
       yy(2) = 0.0D0
       zz(2) = 0.0D0
       ww(2) = ww(1)

       xx(3) = 0.0D0
       yy(3) = 1.0D0
       zz(3) = 0.0D0
       ww(3) = ww(1)

       xx(4) = 0.0D0
       yy(4) =-1.0D0
       zz(4) = 0.0D0
       ww(4) = ww(1)

       xx(5) = 0.0D0
       yy(5) = 0.0D0
       zz(5) = 1.0D0
       ww(5) = ww(1)

       xx(6) = 0.0D0
       yy(6) = 0.0D0
       zz(6) =-1.0D0
       ww(6) = ww(1)

       xx(7) = srti
       yy(7) = srti
       zz(7) = 0.0D0
       ww(7) = 1.0/15.0

       xx(8) =-srti
       yy(8) = srti
       zz(8) = 0.0D0
       ww(8) = ww(7)

       xx(9) = srti
       yy(9) =-srti
       zz(9) = 0.0D0
       ww(9) = ww(7)

       xx(10) =-srti
       yy(10) =-srti
       zz(10) = 0.0D0
       ww(10) = ww(7)

       xx(11) = 0.0D0
       yy(11) = srti
       zz(11) = srti
       ww(11) = ww(7)

       xx(12) = 0.0D0
       yy(12) =-srti
       zz(12) = srti
       ww(12) = ww(7)

       xx(13) = 0.0D0
       yy(13) = srti
       zz(13) =-srti
       ww(13) = ww(7)

       xx(14) = 0.0D0
       yy(14) =-srti
       zz(14) =-srti
       ww(14) = ww(7)

       xx(15) = srti
       yy(15) = 0.0D0
       zz(15) = srti
       ww(15) = ww(7)

       xx(16) =-srti
       yy(16) = 0.0D0
       zz(16) = srti
       ww(16) = ww(7)

       xx(17) = srti
       yy(17) = 0.0D0
       zz(17) =-srti
       ww(17) = ww(7)

       xx(18) =-srti
       yy(18) = 0.0D0
       zz(18) =-srti
       ww(18) = ww(7)

c-----------
      end if
c-----------

      Return
      End
