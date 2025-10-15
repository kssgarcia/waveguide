      subroutine chan2l_cd
     +
     + (N      ! number of interface segments (elements)
     + ,c      ! concentration at nodes
     + ,c0     ! concentraion at mid-nodes
     + ,Ds     ! surfactant diffusivity
     + ,Dt     ! time step
     + ,Move
     + ,X,Y    ! surface nodes
     + ,crv    ! curvature
     + ,Ux,Uy
     + ,Un,Ut
     + ,srfam0 ! initial amount of surfactant
     + ,srfam  ! total amount of surfactant
     + )

c=======================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=======================================

c-----------------------------------
c Update the surfactant concentration
c over a periodic interface
c using an implicit finite-volume method
c
c SYMBOLS:
c -------
c	
c udn ... inner product of u and n 
c udt ... inner product of u and t
c alen ... arc length of an element
c
c amat ... surfactant concentration coefficient matrix	
c rhs ... right hand side vector of amat matrix equation	
c sol ... solution vector from amat matrix equation	
c 	
c cnt ... temporary surfactant concentration at nodes
c cmt ... temporary surfactant concentration at elements
c
c Move = 0 points move with total velocity
c        1 points move with normal velocity
c
c crv:  curvature
c------------------------------------------------------	

      Implicit Double Precision (a-h,o-z)

      Dimension   X(0:513),Y(0:513)

      Dimension  Ux(0:513),Uy(0:513)
      Dimension  Un(0:513),Ut(0:513)
      Dimension crv(0:513)

      Dimension alen(0:900),udn(0:900),udt(0:900)

      Dimension Amat(900,900),rhs(900),sol(900)

      Dimension   c(0:513)
      Dimension  c0(0:513)

      Dimension ct(0:900),c0t(0:900)  ! temporary storage

c--------
c prepare
c--------

      Ds2 = 2.0D0*Ds

      N1 = N+1

c------------------------------
c compute preliminary variables
c at the nodes
c------------------------------

      Do i=1,N
         udn(i) = -Un(i)
         udt(i) =  Ut(i)
        alen(i) =  Dsqrt((X(i+1)-X(i))**2
     +                  +(Y(i+1)-Y(i))**2)
c        write (6,*) udt(i)
      End Do 

c------------
c wrap around
c------------

       udn(0) =  udn(N)
       udt(0) =  udt(N)
      alen(0) = alen(N)

       udn(N1) =  udn(1)
       udt(N1) =  udt(1)
      alen(N1) = alen(1)

c-------------------------------------
c initialize a temporary concentration
c-------------------------------------

      Do i=0,N
       c0t(i) = 0.0D0
      End Do

c----------------------------------
c Generate a finite-volume matrix
c by the method of impulses
c----------------------------------

      Do 1 i=1,N

      c0t(i) = 1.0D0   ! impulse at the ith mid-point

c---
c wrap around
c---

      c0t(0)  = c0t(N)
      c0t(N1) = c0t(1)

c---
c concentration at the surface nodes
c by interpolation
c---

      Do j=1,N
        ja = j-1
        aa = alen(ja)
        bb = alen(j)
        ct(j) = (c0t(ja)*bb+c0t(j)*aa)/(aa+bb)
      End Do

      ct(0)  = ct(N)
      ct(N1) = ct(1)

c---
c Generate the ith column of the matrix
c---

      Do j=1,N

       ja = j-1
       j1 = j+1

       Amat(j,i) = 
     +    - ct(j1)*udt(j1)+ct(j)*udt(j)
     +    - c0t(j)*0.5D0*(udn(j) *crv(j)
     +                   +udn(j1)*crv(j1))*alen(j)
     +    +Ds2*(c0t(j1)-c0t(j) )/(alen(j1)+alen(j))
     +    -Ds2*(c0t(j) -c0t(ja))/(alen(j) +alen(ja))

c-------
c points move with total velocity:
c------

       if(Move.eq.0) then
        Amat(j,i) = Amat(j,i)
     +            + 0.5D0*(ct(j1)-ct(j))*(udt(j1)+udt(j))
       end if

c------

       Amat(j,i) = -Amat(j,i)*Dt/alen(j)

      End Do

      Amat(i,i) = 1.0D0+Amat(i,i)

      c0t(i) = 0.0D0   ! reset

  1   Continue

c-----------------
c Generate the RHS
c-----------------

      Do i=1,N
       rhs(i) = c0(i)
      End do

c-------------------
c call matrix solver
c (gauss elimination)
c-------------------

      Isym_g = 0    ! system is not symmetric
      Iwlpvt = 1    ! pivoting enabled

      call gel
     +
     +   (N    ! system size
     +   ,amat
     +   ,rhs
     +   ,sol
     +   ,Isym_g
     +   ,Iwlpvt
     +   ,det
     +   ,Istop
     +   )

c--------
c display
c--------

      Isee = 0

      if(Isee.eq.1) then

       Do i=1,N
        write (6,109) (Amat(i,j),j=1,N),rhs(i),sol(i)
       End Do

      end if

c-----------------------------
c assign solution to mid-nodes
c-----------------------------

      Do i=1,N
        c0(i) = sol(i)
      End Do 

c-----
c wrap
c-----

      c0(0)  = c0(N)
      c0(N1) = c0(1)

c---------------------------------------
c compute the total amount of surfactant
c---------------------------------------

      srfam = 0.0D0
      altot = 0.0D0

      Do i=1,N
       srfam = srfam + c0(i)*alen(i)
       altot = altot +       alen(i)
      End Do

c     write (6,200) srfam

c---------------------------------
c shift to conserve the surfactant
c---------------------------------

      shift = (srfam0-srfam)/altot
c     shift = 0.0D0

      Do i=1,N
       c0(i) = c0(i) + shift
      End Do

c-----------------------------------
c concentration at the surface nodes
c by interpolation
c----------------------------------

      Do j=1,N
        ja = j-1
        aa = alen(ja)
        bb = alen(j)
        c(j) = (c0(ja)*bb+c0(j)*aa)/(aa+bb)
      End Do

c-----
c wrap
c-----

      c(N1) = c(1)
      c(0)  = c(N)

c     Do i=0,N1
c       write (6,100) i,c(i),Ut(i),Un(i)
c     End Do
c     pause

c-----
c Done
c-----

 200  Format (" Total amount of surfactant ",F15.10)
 100  Format (1X,I3,10(1X,f12.8))
 103  Format (1X,I3,1X,I3,10(1X,f12.8))
 109  Format (500(1X,f5.3))

      return
      end
