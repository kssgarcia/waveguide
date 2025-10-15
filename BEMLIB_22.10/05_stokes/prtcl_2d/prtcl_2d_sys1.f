      subroutine prtcl_2d_sys1
     +
     +   (Iflow
     +   ,expn
     +   ,Niter
     +   ,tol
     +   ,Istop
     +   )

c-------------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c------------------------------------------

c------------------------------------------------
c Compute:
c
c  the particle influence matrices (PIM)
c  the right-hand side of the linear system (PRHS)
c---------------------------------------------

c------------------------------------------
c SYMBOLS:
c --------
c
c PIM(i,j,NE(i),NE(j)) :
c
c    particle interaction matrix for particles
c    i and j
c
c NE(i): number of boundary elements around
c        particle i
c
c pev:  particle eigenvector
c pevt: particle transpose eigenvector
c
c CAPACITY:
c ---------
c
c    72 particles
c    64 elements along each particle
c
c------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  NE(72),Itp(72)

      Dimension axis1(72),axis2(72)
      Dimension xcntr(72),ycntr(72),tilt(72)

      Dimension  xw(72,65),yw(72,65),tw(72,65)
      Dimension  fx(72,64),fy(72,64)

      Dimension   X0(4608),  Y0(4608),T0(4608)
      Dimension  ux0(4608), uy0(4608)
      Dimension elml(4608)
      Dimension vnx0(4608),vny0(4608)
      Dimension tnx0(4608),tny0(4608)

      Dimension link(72,72)

      Dimension pev(128), pevt(128)

      Dimension PIM(49,49,128,128)
      Dimension PRHS(72,128)

c----
c for solving the master linear system (MLS)
c----

      Dimension Nbl(72),Nsz(72)

      Dimension Ifill(72,72)

      Dimension Ablock(500,500),Ainvblock(500,500)
      Dimension Bblock(500),SOLblock(500)

      Dimension SOL(4608),SOLsave(4608)

      Dimension BLOCKinv(72,500,500)

c--------------
c common blocks
c--------------

      common/INTGR1/Nprtcl,NE,Itp,Ncl
      common/INTGR2/NGL

      common/pax/axis1,axis2
      common/pap/xcntr,ycntr,tilt

      common/points/xw,yw,tw

      common/colloc1/ux0,uy0
      common/colloc2/vnx0,vny0,elml
      common/colloc3/tnx0,tny0
      common/colloc4/x0,y0,t0
      common/colloc6/fx,fy

      common/piii/pi,pi2,pi4

c----------
c constants
c----------

      Null = 0

c-----------
c Initialize
c-----------

      Istop = 0

c-------------------------
c Generate the link matrix
c-------------------------

      call Ilink
     +
     +   (Nprtcl
     +   ,Iflow
     +   ,expn
     +   ,link
     +   )

c---
c print
c---

      write (6,*)
      write (6,*) " prtcl_2d_sys1: link matrix:"
      write (6,*)

      Do i=1,Nprtcl
       write (6,*) (link(i,j),j=1,Nprtcl)
      End Do

      pause

c-------------------------------
c Generate: 
c
c   the influence matrices
c   consisting of integrals of the
c   single-layer potential
c
c   the right-hand side
c------------------------------
	
      write (6,*) " prtcl_2d_sys1: Generating the particle"
      write (6,*) "                interaction matrices"

      Iopt = 1     ! need only the slp

      Ic = 0       ! collocation point counter

      Do 31 i=1,Nprtcl     ! loop over particles

       write (6,*) " prtcl_2d_sys1: Collocating over particle ",i

       Do 32 J=1,NE(i)    ! loop over particle collocation points

       Ic = Ic+1

c      write (6,*) " prtcl_2d_sys1: Collocating at point :",Ic

       js = j+NE(i)

c----
c particle eigenvectors
c----

       pev(j)  = vnx0(Ic)
       pev(js) = vny0(Ic)

       pevt(j)  = vnx0(Ic)*elml(Ic)
       pevt(js) = vny0(Ic)*elml(Ic)

c---
c right-hand side
c---

       PRHS(i,j ) = -pi4*ux0(Ic) 
       PRHS(i,js) = -pi4*uy0(Ic) 

c      sum5 = 0.0D0   ! for debugging
c      sum6 = 0.0D0

       Kc = 0       ! element counter

c---
c particle influence matrices
c---

       Do 33 k=1,Nprtcl       ! loop over particles
        Do 34 l=1,NE(k)       ! loop over elements

         kc = kc+1

         X1 = XW(k,l)          ! first element point
         Y1 = YW(k,l)
         T1 = TW(k,l)

         X2 = XW(k,l+1)        ! second element point
         Y2 = YW(k,l+1)
         T2 = TW(k,l+1)

         Ising = 0
         If(Ic.eq.kc) Ising = 1

         call prtcl_2d_sdlp
     +
     +     (Iopt
     +     ,Iflow
     +     ,X0(Ic),Y0(Ic),t0(Ic)
     +     ,X1,Y1,T1
     +     ,X2,Y2,T2
     +     ,NGL
     +     ,Ising
     +     ,Itp(k)
     +     ,axis1(k),axis2(k)
     +     ,xcntr(k),ycntr(k),tilt(k)
     +     ,QQxx,QQxy
     +     ,QQyx,QQyy
     +     ,Wxx,Wyx
     +     ,Wxy,Wyy
     +     )

         ls = l+NE(k)

         PIM(i,k,j, l)  = QQxx     ! interaction of i and k particle
         PIM(i,k,j, ls) = QQyx
         PIM(i,k,js,l)  = QQxy
         PIM(i,k,js,ls) = QQyy

c        sum5 = sum5 + QQxx*vnx0(Kc)+ QQyx*vny0(Kc)
c        sum6 = sum6 + QQxy*vnx0(Kc)+ QQyy*vny0(Kc)

   34     Continue   ! end of loop over elements
   33    Continue   ! end of loop over particles

c---------
c  It should be: sum5 = 0, sum6 = 0
c---
c     write (6,112) sum5,sum6
c---------

  32   Continue    ! end of loop over particle collocation points

c--------------------------------
c Regularization by deflation
c deflate the i-th particle 
c self-interaction matrix
c--------------------------------

       Do j=1,NE(i)

        js = j+NE(I)

        Do l=1,NE(i)
         ls = l+NE(I)
           PIM(i,i,j, l)  = PIM(i,i,j, l)  + pev(j) *pevt(l)
           PIM(i,i,j, ls) = PIM(i,i,j, ls) + pev(j) *pevt(ls)
           PIM(i,i,js,l)  = PIM(i,i,js,l)  + pev(js)*pevt(l)
           PIM(i,i,js,ls) = PIM(i,i,js,ls) + pev(js)*pevt(ls)
        End Do

       End Do

  31  Continue    ! end of loop over particles

c------------------------------
c Done generating the 
c particle interaction matrices
c------------------------------

c-----------------------------
c Generate and compute the inverse matrices
c of the particle clouds
c-----------------------------

      write (6,*) " prtcl_2d_sys1:"
      write (6,*) 
      write (6,*) "     Defining and computing"
      write (6,*) "     the inverse matrices"
      write (6,*) "     of the particle cloud matrices"

      Isym_g = 0   ! system is not symmetric
      Iwlpvt = 0   ! pivoting disabled

c-----------
      Do I=1,Nprtcl    ! loop over particle clouds
c-----------

c----
c identify the ith cluster
c
c Nsz(i):
c
c    number of particles in the
c    ith particle-cluster
c
c Ifill(i,j):
c
c    label of jth particle in the cluster,
c    where j=1,...,Nsz(i)
c
c Nbl(i): size of the ith particle-cluster
c----

      Nsz(i)     = 1
      Ifill(i,1) = I
      Nbl(i)     = 2*NE(i)

      Do K=1,Nprtcl
        If(link(i,K).eq.1) then

         Nsz(i) = Nsz(i)+1
         Ifill(i,Nsz(i)) = K
         Nbl(i) = Nbl(i)+2*NE(K)

        End If
      End Do

c----
c generate the ith particle cloud matrix:
c
c Ablock
c----

      Im = 0                 ! row counter

      Do Icl=1,Nsz(I)        ! run over particle clusters

        nl = Ifill(i,Icl)    ! particle label

        Do J=1,2*NE(nl)

         Im = Im+1
         Jm = 0                 ! column counter

          Do Jcl=1,Nsz(I)       ! run over particles

           ml = Ifill(I,Jcl)

            Do L=1,2*NE(ml)
             Jm = Jm+1
             Ablock(Im,Jm) = PIM(nl,ml,J,L)
            End Do

          End Do

        End Do

      End Do

c---
c compute the inverse
c---

      write (6,*) " prtcl_2d_sys1: computing the inverse of cloud", i
      write (6,*) "                matrix size", Nbl(i)

      call gel_inv    ! compute the inverse
     +
     +   (Nbl(i)
     +   ,Ablock
     +   ,Ainvblock
     +   ,Isym_g
     +   ,Iwlpvt
     +   ,Deter
     +   ,Istop
     +   )

      Do J=1,Nbl(I)
        Do L=1,Nbl(I)
           BLOCKinv(I,J,L) = Ainvblock(J,L)
        End Do
      End Do

c-----------
      End Do    ! end of loop over particles
c-----------

c------------------------------------
c Done computing the inverse matrices
c------------------------------------

c---------------------
c Start the iterations
c---------------------

      Do Iter=1,Niter

c--- PRINTING
      Go to 881
      Ic = 0 
      Do i=1,Nprtcl
        write (6,*) NE(i),i
        Do J=1,NE(i)
         Ic = Ic+1
         write (6,102) J,X0(Ic),Y0(Ic),fx(I,J),fy(I,J)
     +                  
        End Do
      End Do
 881  Continue
c---

c---
c fill in the SOL vector
c---

      Ic = 0
      Do I=1,Nprtcl
       Do J=1,NE(I)
        Ic = Ic+1
        SOL(Ic) = fx(I,J)
       End Do
       Do J=1,NE(I)
        Ic = Ic+1
        SOL(Ic) = fy(I,J)
       End Do
      End Do

c-----------
      Do I=1,Nprtcl    ! loop over particle clouds

c---
c compute the RHS
c---

      Im = 0       !  point counter

      Do Icl=1,Nsz(I)

        nl = Ifill(I,Icl)      ! particle label

        Do J=1,2*NE(nl)

         Im = Im+1

         Bblock(Im) = PRHS(nl,J) 

          Do k=1,Nprtcl

           If(link(I,k).eq.0) then
            Do l=1,NE(k)
             ls = l+NE(k)
              Bblock(Im) = Bblock(Im) - PIM(nl,k,J,l) *fx(k,l)
     +                                - PIM(nl,k,J,ls)*fy(k,l)
            End Do
           End If

          End Do

        End Do
      End Do

c---
c projection in terms of the inverse
c---

       Do J=1,Nbl(I)

        SOLblock(J) = 0.0D0

        Do K=1,Nbl(I)
         SOLblock(J) = SOLblock(J)+BLOCKinv(I,J,K)*Bblock(K)
        End Do

       End Do

c---------------------
c update the tractions
c by dismantling the traction blocks
c---------------------

      Im = 0   ! element counter

      Do Icl=1,Nsz(I)

       nl = Ifill(I,Icl)

       Do J=1,NE(nl)
         JS = J+NE(nl)
         fx(nl,J) = SOLBlock(Im+J)
         fy(nl,J) = SOLBlock(Im+JS)
       End Do

       Im = Im+2*NE(nl)

      End Do

c-----------
      End Do    ! end of loop over clouds
c-----------

c--------------------------
c compute the global solution vector
c and the error
c---------------------------

      Ic = 0

      Do I=1,Nprtcl

       Do J=1,NE(I)
        Ic = Ic+1
         SOLsave(Ic) = Sol(Ic)
         SOL(Ic) = fx(I,J)
       End Do

       Do J=1,NE(I)
        Ic = Ic+1
         SOLsave(Ic) = Sol(Ic)
         SOL(Ic) = fy(I,J)
       End Do

      End Do

      Error  = 0.0D0
      Errmax = 0.0D0

       Do j=1,2*Ncl
        error = error+(SOLsave(j)-SOL(j))**2
        Diff = abs(SOLsave(j)-SOL(j))
        If(Diff.gt.Errmax) Errmax = Diff
       End Do

      error = sqrt(error)/(2*Ncl)

      write (6,888) Iter,Errmax

      If(error.le.tol) Go to 99

c------------
       End Do ! loop over iterations
c------------

 99   Continue

c-----
c Done
c-----

 101  Format (80(1x,f7.3))
 102  Format (1x,i3,20(1x,f10.5))
 103  Format (1x,i3,5(1x,f15.10))
 104  Format (1x,i3,20(1x,f9.5))
 111  Format (80(10(1x,f7.3),/))
 112  Format (80(1x,f10.7))

 888  Format (1x,"Iter: ",i3," Point corr: ",f15.10)

      Return
      End
