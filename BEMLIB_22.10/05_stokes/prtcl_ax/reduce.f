      subroutine reduce ()

c-----------------------------------------
c FDLIB and  BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c--------------------------------------------
c Reduce the primary or preconditioned system
c
c Replace the last equation of each particle block
c with the constraint:
c
c  Sum [ f.n elar ] = 1.0
c
c LEGEND:
c -------
c
c AL:   Master Matrix
c BL:   Master rhs
c--------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension NE(25),Itp(25)
      Dimension vnx0(3200),vny0(3200),elar(3200)
 
      Dimension AL(3200,3200),BL(3200)

c--------------
c common blocks
c--------------

      common/INTGR1/Nprtcl,NGL,NE,Itp,Ncl
      common/colloc2/vnx0,vny0,elar
      common/matrix/AL,BL

c--------
c prepare
c--------

      write (6,*) " reduce: Reducing the MLS"

      Ncl2 = 2*Ncl

c-----------
c initialize
c-----------

      Ir  = 0     ! row pointer
      Ic  = 0     ! column counter
      Icp = 0     ! collocation point counter

c-------
c reduce
c-------

      Do i=1,Nprtcl

       Ir = Ir+2*NE(i)    ! last row of ith particle block

       Do j=1,Ncl2
        AL(Ir,j) = 0.0D0   ! zero out
       End Do

       BL(Ir) = 1.0D0    ! rhs

       Do j=1,NE(i)
        Icp = Icp+1
        Icj = Ic+j
        AL(Ir,Icj) = vnx0(Icp)*elar(Icp)
        AL(Ir,Icj) = vnx0(Icp)
        Icj = Icj+NE(i)
        AL(Ir,Icj) = vny0(Icp)*elar(Icp)
        AL(Ir,Icj) = vny0(Icp)
       End Do

       Ic = Ic+2*NE(i)

      End Do

c-----
c Done
c-----

      write (6,*) " reduce: done"

 112  Format (80(1x,f10.7))

      Return
      End
