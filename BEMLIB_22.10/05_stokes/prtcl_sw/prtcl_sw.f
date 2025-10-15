      program prtcl_sw

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c---------------------------------------------
c Swirling flow due to the axial rotation of an 
c axisymmetric particle.
c
c This program solves
c an integral equation for the boundary traction
c and computes the torque exerted on the particle.
c
c The particle contour in the xy plane,
c consists of a Number of SeGments (NSG)
c
c The segments may be straight lines or circular arcs.
c
c Each segment is discretized into a specified
c number of elements.
c
c
c Symbols:
c --------
c
c NSG: Number of segments defining the particle contour
c
c NE(i): Number of elements on ith segment 
c
c RT(i): Stretch ratio of elements on ith segment 
c
c Itp(i): Index for shape of the ith segment:
c         1 for a straight segment
c         2 for a circular arc
c
c (Xe, Ye):  end-nodes of elements on a segment
c (Xm, Ym):  mid-nodes of elements on a segment
c (Xw, Yw):  end-nodes of elements on all segments
c
c ff:	phi-component of the traction at elements
c
c (X0, Y0):  coordinates of collocation points
c
c t0(i): angle subtended from the center of a circular element
c        at ith collocation point
c
c arel(i):  axisymmetric surface area of ith element
c
c NGL: Number of Gaussian points for integration 
c      over each element
c
c
c Note:
c -----
c
c The normal vector points into the flow
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension NE(10),Itp(10)
      Dimension actis(10),xcntr(10),ycntr(10)

      Dimension xw(10,200),yw(10,200),tw(10,200)
      Dimension ff(10,200)

      Dimension X0(500),Y0 (500),t0(500),s0(0:500)
      Dimension arel(500)

      Dimension AL(500,500),BL(500),SOL(500)  ! for the linear system

      Dimension ZZ(20),WW(20)

c--------------
c common blocks
c--------------

      common/piii/pi,pih,pi2,pi4

      common/XF1/visc,Omega
      common/XF2/NSG,NE,Itp,NGL
      common/XF3/xw,yw,tw
      common/XF4/X0,Y0,T0,S0
      common/XF5/arel
      common/XF5/actis,xcntr,ycntr

      common/ZZWW/ZZ,WW

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pih = 0.5D0*pi
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi

      Null = 0
      None = 1
      Ntwo = 2

      zero  = 0.0D0
      five  = 5.0D0
      fivem =-5.0D0

c-----------------
c Return to repeat
c-----------------

 98   Continue

      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 0 to quit"
      write (6,*) " 1 for a sphere"
      write (6,*) " 2 for a spheroid"
      write (6,*) " 3 for a triangular torus"
      write (6,*) " ------------------------"

      read (5,*) Iflow

      If(Iflow.eq.0) Go to 99

c----------------------
c  Element distribution
c----------------------

      call prtcl_sw_geo (Iflow,Ncl)

c-----------------------
c Printing session
c
c print particle contour 
c-----------------------

      open (9,file="prtcl_sw.xy")

      Do j=1,NSG

       write (9,102) NE(j)+1
       write (6,102) NE(j)+1

       Do i=1,NE(j)+1
        write (9,102) i,XW(j,i),YW(j,i)
        write (6,102) i,XW(j,i),YW(j,i)
       End Do

      End Do

c---
c Display collocation points
c---

      write (6,*)
      write (6,*) Ncl," Collocation points "
      write (6,*)

      write (6,*)
      write (6,*) "Collocation points"
      write (6,*)

      Do i = 1,Ncl
        write (6,102) i,X0(i),Y0(i)
      End Do

c---------------
c  PREPARE TO RUN
c----------------

      call gauss_legendre (NGL,ZZ,WW)

c----------------------------------------------
c  Generate the linear system for the traction
c  at the collocation points
c
c  The influence matrix consists
c  of integrals of the single-layer potential
c
c  Compute the rhs by evaluating the dlp
c---------------------------------------------
	
      Do 31 I=1,Ncl    ! loop over collocation points

c      write (6,*)        " Collocation point :",I

       BL(I) = - visc*Omega*Y0(I)     ! right-hand side

       J = 0              ! counter

       Do k=1,NSG         ! loop over segments

        If(Itp(k).eq.2) then
         rad  = actis(k)
         xcnt = xcntr(k)
         ycnt = ycntr(k)
        End If

        Do L=1,NE(k)       ! loop over elements

         X1 = XW(K,L)
         Y1 = YW(K,L)
         T1 = TW(K,L)

         X2 = XW(K,L+1)
         Y2 = YW(K,L+1)
         T2 = TW(K,L+1)

         J = J+1
         Ising = 0
         If(I.eq.J) Ising = 1

         call prtcl_sw_slp
     +
     +    (X0(I),Y0(I),t0(I)
     +    ,X1,Y1,T1
     +    ,X2,Y2,T2
     +    ,NGL
     +    ,Ising
     +    ,Itp(k)
     +    ,rad,xcnt,ycnt
     +    ,QQQ
     +    )

         AL(I,J) = QQQ

         End Do

        End Do

  31  Continue

c------------------------
c SOLVE THE LINEAR SYSTEM
c------------------------

      write (6,*)
      write (6,*) "Size of the linear system:",Ncl
      write (6,*)

      Isym_g = 0    ! system is not symmetric
      Iwlpvt = 1    ! pivoting enabled

      call gel
     +
     +   (Ncl
     +   ,AL,BL,SOL
     +   ,Isym_g
     +   ,Iwlpvt
     +   ,Deter
     +   ,Istop
     +   )

c----------------------
c DISPLAY LINEAR SYSTEM
c----------------------

c     Do I=1,Ncl
c      write (6,101) (AL(I,J),J=1,Ncl),BL(I),SOL(I)
c     End Do

c---------------------
c  DISTRIBUTE SOLUTION
c---------------------

      k = 0        ! counter

      Do i=1,NSG
       Do j=1,NE(i)
        k = k+1
        ff(i,j) = SOL(k)
       End Do
      End Do

c---------------
c print solution
c---------------

      open (3,file="prtcl_sw.trc")

      k = 0

      write (6,*)
      write (6,*)  " x, y, arc length, traction"
      write (6,*)

      Do i=1,NSG

        write (6,102) NE(i)
        write (3,102) NE(i)

        Do j=1,NE(i)
         k = k+1
         write (6,102) J,X0(k),Y0(k),s0(k),ff(i,j)
         write (3,102) J,X0(k),Y0(k),s0(k),ff(i,j)
        End Do

      End Do

c-------------------
c Compute the torque
c-------------------

      torquex = 0.0D0

      k = 0  ! element counter

      Do i=1,NSG
        Do j=1,NE(i)
          k = k+1
          torquex = torquex + Y0(k)*arel(k)*ff(i,j)
        End Do
      End Do

      write (6,*)
      write (6,888) torquex

c-------
c Repeat
c-------

      Go to 98

  99  Continue

c----------
c Finish up
c----------

      write (9,104) Null
      close (9)
      write (3,102) Null
      close (3)

c-----
c Done
c-----

 101  Format (20(1x,f7.3))
 102  Format (1x,i3,20(1x,f15.10))
 104  Format (1x,i3,20(1x,f9.5))
 150  Format (1X,10(1X,f10.5))
 888  Format (1X,"Axial Torque: ",f10.5)

      Stop
      End
