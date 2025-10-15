      subroutine caps_2d_cd
     +
     + (Move
     + ,Isym
     + ,vnx,vny
     + ,vtx,vty
     + ,cm
     + ,c
     + ,crv
     + ,Ds,Dt
     + ,srfam
     + ,Istop
     + )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c-----------------------------------
c update the surfactant concentration
c over the interface
c by an explicit  finite-volume method
c
c SYMBOLS:
c -------
c	
c udn .... inner product of u and n 
c udt .... inner product of u and t
c alen .... arc length of an element
c
c amat ... surfactant concentration coefficient matrix	
c rhs .... right hand side vector of amat matrix equation	
c xs ..... solution vector from amat matrix equation	
c 	
c u,v ..... x,y velocity components of material surface	 
c vnx,vny .. x and y components of surface normal vector	
c vtx,vty .. x and y components of surface tangent vector	
c s .... arc length	
c c .... surfactant concentration at nodes
c cm .... surfactant concentration at elements
c Ds ... surfactant diffusivity	
c nsg ... number of surface segments
c	
c cnt .... temporary surfactant concentration at nodes
c cmt .... temporary surfactant concentration at elements
c
c Move = 0 points move with total velocity
c        1 points move with normal velocity
c
c Ijac... size of the Jacobian
c------------------------------------------------------	

      Implicit Double Precision (a-h,o-z)

      Dimension   u(0:900),  v(0:900)
      Dimension   c(0:900), cm(0:900)
      Dimension cmt(0:900),cnt(0:900)  ! remporary storage

      Dimension  XC(900), YC(900),  R(900),   S(900)
      Dimension th1(900),th2(900),th3(900),ornt(900)

      Dimension vnx(900),vny(900),vtx(900),vty(900)
      Dimension crv(900)

      Dimension alen(0:900),udn(0:900),udt(0:900)

      Dimension amat(500,500),rhs(500),xs(500)

c--------------
c common blocks
c--------------

      common/UUVV/u,v
      common/ARCC/XC,YC,R,S,TH1,TH2,TH3,ORNT

      common/ancI1/NSG,NSG1,NSG2,NSGM,NSGQ

      common/piii/pi,pih,piq,pi2,pi4,pi6,pi8,srpi

c----------
c constants
c----------

      pi = 3.14159 265358 979 32384D0

c--------
c prepare
c--------

      Ds2 = 2.0D0*Ds

      if(Isym.eq.0) then
        Ijac = NSG
      else if(Isym.eq.2) then
        NSGQ = NSG/4+1
        NSGM = NSG/2+1
        Ijac = NSGQ-1
      end if
     
c------------------------------
c compute preliminary variables
c------------------------------

      Do i=1,NSG
         udn(i) = u(i)*vnx(i)+v(i)*vny(i)  
         udt(i) = u(i)*vtx(i)+v(i)*vty(i)
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
c generate the finite-volume matrix
c by the method of impulses
c----------------------------------

      Do i=1,Ijac

      cmt(i) = 1.0   ! impulse at the ith point

c-----------------------
      if(Isym.eq.2) then
        Do k=1,NSGQ
          k1      = NSGM-k+1
          cmt(k1) =  cmt(k)
          k1      = NSGM+k-1
          cmt(k1) =  cmt(k)
          k1      = NSG1-k+1
          cmt(k1) =  cmt(k)
        End Do
      end if
c-----------------------

c---
c wrap around
c---

      cmt(0)    = cmt(NSG)
      cmt(NSG1) = cmt(1)

c---
c concentration at end-points
c by interpolation
c---

      Do j=1,NSG
        ja = j-1
        aa = alen(ja)
        bb = alen(j)
        cnt(j) = (cmt(ja)*bb+cmt(j)*aa)/(aa+bb)
      End Do

c---
c wrap around
c---

      cnt(0)    = cnt(NSG)
      cnt(NSG1) = cnt(1)

c---
c generate the ith column of the matrix
c---

      Do j=1,Ijac

       ja = j-1
       j1 = j+1

       amat(j,i) = 
     +    - cnt(j1)*udt(j1)+cnt(j)*udt(j)
     +    - cmt(j)*0.5D0*(udn(j)*crv(j)+udn(j1)*crv(j1))*alen(j)
     +    +Ds2*(cmt(j1)-cmt(j) )/(alen(j1)+alen(j))
     +    -Ds2*(cmt(j) -cmt(ja))/(alen(j) +alen(ja))

       if(Move.eq.0) then
       amat(j,i) = amat(j,i)
     +           + 0.5D0*(cnt(j1)-cnt(j))*(udt(j1)+udt(j))
       end if

       amat(j,i) = -amat(j,i)*dt/alen(j)

      End Do

      amat(i,i) = 1.0D0+amat(i,i)

      cmt(i) = 0.0D0 ! reset

      End Do

      Do i=1,Ijac
       rhs(i)= cm(i)
      End Do

c-------------------
c call matrix solver
c-------------------

      Isym_g = 0    ! system is not symmetric
      Iwlpvt = 1    ! pivoting enabled

      call gel
     +
     +  (Ijac  ! system size
     +  ,amat
     +  ,rhs,xs
     +  ,Isym_g,Iwlpvt
     +  ,det
     +  ,Istop
     +  )

c-----------------------------
c assign solution to mid-nodes
c-----------------------------

      Do i=1,Ijac
        cm(i) = xs(i)
      End Do 

c------------------------------
c reflect a symmetric interface
c------------------------------

      if(Isym.eq.2) then
        Do i=1,NSGQ
          i1     = NSGM-i+1
          cm(i1) =  cm(i)
          i1     = NSGM+i-1
          cm(i1) =  cm(i)
          i1     = NSG1-i+1
          cm(i1) =  cm(i)
        End Do
      end if

c-----
c wrap
c-----

      cm(0)     = cm(NSG)
      cm(NSG+1) = cm(1)
      cm(NSG+2) = cm(2)
      cm(NSG+3) = cm(3)
      cm(NSG+4) = cm(4)
      cm(NSG+5) = cm(5)

c-----------------------------------
c compute total amount of surfactant
c-----------------------------------

      srfam = 0.0D0

      Do i=1,nsg
       srfam = srfam+cm(i)*alen(i)
      End Do

      write (6,200) srfam

c-----
c done
c-----

 200  Format (" Total amount of surfactant ",F15.10)

      Return
      End
