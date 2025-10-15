      subroutine em_2d_cd
     +
     +  (Idrop
     +  ,Move
     +  ,NSG
     +  ,vnx,vny
     +  ,vtx,vty
     +  ,cm
     +  ,c
     +  ,crv
     +  ,Ds,Dt
     +  ,srfam
     +  ,Istop
     +  )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c--------------------------------------------
c update the surfactant concentration
c over an interface
c by an implicit finite-volume method
c
c SYMBOLS:
c -------
c
c Ds ... surfactant diffusivity	
c nsg ... number of interface segments
c
c u,v .....  velocity components at the surface	 
c vnx,vny .. x and y components of surface normal vector	
c vtx,vty .. x and y components of surface tangent vector	
c crv...     curvature at end-nodes
c
c s .... arc length	
c c .... surfactant concentration at end-nodes
c cm .... surfactant concentration at mid-nodes
c
c Move = 0 points move with total velocity
c        1 points move with normal velocity
c
c srfam.... total amount of a surfactant
c
c udn .... inner product of u and n 
c udt .... inner product of u and t
c alen .... arc length of an element
c
c cnt .... temporary surfactant concentration at end-nodes
c cmt .... temporary surfactant concentration at mid-nodes
c
c amat ... surfactant concentration coefficient matrix	
c rhs .... right hand side vector of amat matrix equation	
c sol ..... solution vector from amat matrix equation	
c------------------------------------------------------	

      Implicit Double Precision (a-h,o-z)

      Dimension   u(0:200),  v(0:200)
      Dimension   c(0:200), cm(0:200)
      Dimension cmt(0:200),cnt(0:200)    ! temporary storage

      Dimension  XC(200), YC(200),  R(200),   S(200)
      Dimension TH1(200),TH2(200),TH3(200),ORNT(200)
 
      Dimension vnx(0:200),vny(0:200)
      Dimension vtx(0:200),vty(0:200)
      Dimension crv(0:200)

      Dimension alen(0:200),udn(0:200),udt(0:200)

      Dimension amat(200,200),rhs(200),sol(200)

c--------------
c common blocks
c--------------

      common/UUVV/u,v
      common/ARCC/XC,YC,R,S,TH1,TH2,TH3,ORNT
      common/piii/pi,pih,piq,pi2,pi4,pi6,pi8,srpi

c--------
c prepare
c--------

      NSG1 = NSG+1

      Ds2 = 2.0D0*Ds

c------------------------------
c compute preliminary variables
c------------------------------

      Do i=1,NSG
         udn(i) = u(i)*vnx(i) + v(i)*vny(i)  
         udt(i) = u(i)*vtx(i) + v(i)*vty(i)
        alen(i) = s(i+1)-s(i)
      End Do 

c------------
c wrap around
c------------

       udn(0) =  udn(NSG)
       udt(0) =  udt(NSG)
      alen(0) = alen(NSG)

       udn(NSG1) =  udn(1)
       udt(NSG1) =  udt(1)
      alen(NSG1) = alen(1)

c-----------
c initialize
c-----------

      Do i=0,NSG1
       cmt(i) = 0.0D0
      End Do

c----------------------------------
c generate a finite-volume matrix
c by the method of impulses
c----------------------------------

      Do i=1,NSG

      cmt(i) = 1.0D0  ! impulse

      cmt(0)    = cmt(NSG)  ! wrap
      cmt(NSG1) = cmt(1)    ! wrap

c---
c concentration at end-nodes
c by linear interpolation
c---

      Do j=1,NSG
        ja = j-1
        aa = alen(ja)
        bb = alen(j)
        cnt(j) = (cmt(ja)*bb+cmt(j)*aa)/(aa+bb)
      End Do

c------------
c wrap around
c------------

      cnt(0)    = cnt(NSG)
      cnt(NSG1) = cnt(1)

c---
c Generate the ith column of the matrix
c---

      Do j=1,NSG

       ja = j-1
       j1 = j+1

       amat(j,i) = 
     +    - cnt(j1)*udt(j1)+cnt(j)*udt(j)
     +    - cmt(j)*0.5D0*(udn(j)*crv(j)+udn(j1)*crv(j1))*alen(j)
     +    +Ds2*(cmt(j1)-cmt(j) )/(alen(j1)+alen(j))
     +    -Ds2*(cmt(j) -cmt(ja))/(alen(j) +alen(ja))

       if(move.eq.0) then
       amat(j,i) = amat(j,i)
     +           + 0.5D0*(cnt(j1)-cnt(j))*(udt(j1)+udt(j))
       end if

       amat(j,i) = -amat(j,i)*Dt/alen(j)

      End Do

      amat(i,i) = 1.0D0+amat(i,i)

      cmt(i) = 0.0D0   ! reset

      End Do

      Do i=1,NSG
       rhs(i)=cm(i)
      End Do
 
c--------------------
c  call matrix solver
c--------------------

      Isym_g = 0      ! system is not symmetric
      Iwlpvt = 1      ! will pivot

      call gauss_small 
     +
     +  (NSG
     +  ,amat
     +  ,rhs,sol
     +  ,Isym_g
     +  ,Iwlpvt
     +  ,det
     +  ,Istop
     +  )

c-----------------------------
c assign solution to mid-nodes
c-----------------------------

      Do i=1,NSG
        cm(i) = sol(i)
      End Do 

c------------
c wrap around
c------------

      cm(0)     = cm(nsg)
      cm(nsg+1) = cm(1)
      cm(nsg+2) = cm(2)
      cm(nsg+3) = cm(3)
      cm(nsg+4) = cm(4)
      cm(nsg+5) = cm(5)
  
c-----------------------------------
c compute total amount of surfactant
c-----------------------------------

      srfam = 0.0D0

      Do i=1,nsg
       srfam = srfam+cm(i)*alen(i)
      End Do

      write (6,200) Idrop,srfam

c-----
c done
c-----

 200  Format ("Surfactant: Interface ",I3
     +       ,"Total amount of surfactant: ",F15.10)

      return
      end
