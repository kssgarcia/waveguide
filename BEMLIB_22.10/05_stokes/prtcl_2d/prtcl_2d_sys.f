      subroutine prtcl_2d_sys
     +
     +   (Iflow
     +   ,Method
     +   ,Iprec
     +   ,Ireg
     +   ,Nblocks
     +   ,Lump
     +   ,Niter
     +   ,tol
     +   ,Istop
     +   )

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c---------------------------------
c Generate and solve the extended
c Master Linear system (MLS)
c---------------------------------

c------------------------------------------
c SYMBOLS:
c --------
c
c Imn(i): Low  index of ith particle in master matrix
c Imx(i): High index of ith particle in master matrix
c Isz(i):       size of ith particle in master matrix
c
c
c Capacity:
c ---------
c
c    72 particles
c    64 elements along each particle
c
c------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    NE(72),  Itp(72)
      Dimension axis1(72),axis2(72)
      Dimension xcntr(72),ycntr(72),tilt(72)

      Dimension  xw(72,65),yw(72,65),tw(72,65)
      Dimension  fx(72,64),fy(72,64)

      Dimension   X0(4608),  Y0(4608),T0(4608)
      Dimension  ux0(4608), uy0(4608)
      Dimension elml(4608)
      Dimension vnX0(4608),vnY0(4608)
      Dimension tnX0(4608),tnY0(4608)

c----
c global eigenvectors
c-----

      Dimension ev(4608),evt(4608)

c---
c related to the linear solver
c---

      Dimension Imn(72),Imx(72),Isz(72)
      Dimension Nmn(72),Nmx(72),Nsz(72),Lump(72)

c     Dimension AL(4608,4608),BL(4608),SOL(4608) 

      Dimension AL(2304,2304),BL(2304),SOL(2304) 

      Dimension Ablock(500,500),Ainvblock(500,500)
      Dimension Bblock(500),SOLblock(500) 

      Dimension blockinv(72,128,128) 

c---
c related to iterations (888 maximum)
c---

      Dimension SOLsave(4608,888)

      Dimension error(888)    ! rms error
      Dimension perr(888)     ! pointwise error

c     Dimension SOLext(1600,200)

c--------------
c common blocks
c--------------

      common/INTGR1/Nprtcl,NE,Itp,Ncl
      common/INTGR2/NGL

      common/REAL1/visc,shrt,delta,Uinf,wall
      common/REAL2/Uprtcl,Vprtcl,Aprtcl

      common/CHANNELI/IQPD,NGFWW
      common/CHANNELR/U1,U2,RL,h

      common/pax/axis1,axis2
      common/pap/xcntr,ycntr,tilt

      common/points/xw,yw,tw

      common/colloc1/ux0,uy0
      common/colloc2/vnx0,vny0,elml
      common/colloc3/tnx0,tny0
      common/colloc4/x0,y0,t0
      common/colloc6/fx,fy

      common/matrix/AL,BL

      common/piii/pi,pi2,pi4

c----------
c constants
c----------

      Null = 0

c-----------
c Initialize
c-----------

      Istop = 0

c--------
c prepare
c--------

      Ncl2 = 2*Ncl

c-------------------------------
c Generate the influence matrix
c consisting of integrals of the
c single-layer potential
c
c If Method = 2, compute the right-hand side
c by evaluating the double-layer potential
c------------------------------
	
      Ic = 0       ! collocation point counter
      Ib = 0       ! block counter

      Do 31 i=1,Nprtcl     ! loop over particles
       Do 32 j=1,NE(i)     ! loop over collocation points

       Ic = Ic+1

       write (6,*) " prtcl_2d_sys: Collocating at point :",Ic

       Ibj  = Ib  + j      ! row entry for x component of BIE
       Ibji = Ibj + NE(i)  ! row entry for y component of BIE

c---
c right-hand side
c---

       If(method.eq.1) then

        BL(Ibj ) = -pi4*ux0(Ic) 
        BL(Ibji) = -pi4*uy0(Ic)

       Else If(method.eq.2) then

        BL(Ibj)  = -pi2*ux0(Ic)
        BL(Ibji) = -pi2*uy0(Ic)

       End If

c----------------
c      sum1 = 0.0   ! for debugging
c      sum2 = 0.0
c      sum3 = 0.0
c      sum4 = 0.0
c      sum5 = 0.0
c      sum6 = 0.0
c----------------

       kc = 0        ! element counter
       kb = 0        ! block counter

       Do 33 k=1,Nprtcl       ! loop over particles
        Do 34 l=1,NE(k)       ! loop over elements

         kc = kc+1

         kbl  = kb  + l         ! column entry for x traction
         kblk = kbl + ne(k)     ! column entry for y traction

         X1 = XW(k,L)          ! first element point
         Y1 = YW(k,L)
         T1 = TW(k,L)

         X2 = XW(k,L+1)        ! second element point
         Y2 = YW(k,L+1)
         T2 = TW(k,L+1)

         Ising = 0
         If(Ic.eq.kc) Ising = 1

         call prtcl_2d_sdlp
     +
     +     (Method
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

         AL(Ibj,  kbl)  = QQxx
         AL(Ibj,  kblk) = QQyx
         AL(Ibji, kbl)  = QQxy
         AL(Ibji, kblk) = QQyy

         If(Method.eq.2) then
           BL(Ibj)  = BL(Ibj)  + ux0(kc)*Wxx + uy0(kc)*Wyx
           BL(Ibji) = BL(Ibji) + ux0(kc)*Wxy + uy0(kc)*Wyy
c          sum1 = sum1 + Wxx
c          sum2 = sum2 + Wyx
c          sum3 = sum3 + Wxy
c          sum4 = sum4 + Wyy
         End If

c        sum5 = sum5 + QQxx*vnx0(kc)+ QQyx*vny0(kc)
c        sum6 = sum6 + QQxy*vnx0(kc)+ QQyy*vny0(kc)

   34    Continue

         kb = kb + 2*NE(k)     ! shift by one particle block

   33    Continue

c---------
c
c     sum1 = sum1/pi2
c     sum2 = sum2/pi2
c     sum3 = sum3/pi2
c     sum4 = sum4/pi2
c
c  It should be: sum5 = 0, sum6 = 0
c
c  write (6,112) sum1,sum2,sum3,sum4,sum5,sum6
c
c---------

  32  Continue

      Ib = Ib+2*NE(i)     ! shift by one particle block

  31  Continue

c----------------------------------
c Done generating the linear system
c----------------------------------

      write (6,*) " prtcl_2d_sys: Done with the linear system"

c----------------------------
c precondition and regularize
c----------------------------

      If(Iprec.eq.1) then
         write (6,*) " prtcl_2d_sys: Preconditioning"
         call precondition ()
         write (6,*) " prtcl_2d_sys: Done"
c        pause
      End If

      If(Ireg.eq.1) then
         write (6,*) " prtcl_2d_sys: Regularizing"
         call reduce ()
         write (6,*) " prtcl_2d_sys: Done"
c        pause
      End If

c----------------------------
c regularization by deflation
c----------------------------

c---
      If(Ireg.eq.3) then
c---

       write (6,*) " prtcl_2d_sys: Deflating"

        call eigenvectors (ev, evt)

        Ib = 0

        Do i=1,Nprtcl  ! run over particles

         Min = Ib+1          ! lower limit of ith partition
         Max = Ib+2*NE(i)    ! upper limit of ith partition

         Do j=Min,Max
          Do k=Min,Max
           AL(j,k) = AL(j,k) + ev(j)*evt(k)
          End Do
         End Do

         Ib = Ib+2*NE(i)

       End Do

       write (6,*) " prtcl_2d_sys: Done"
c      pause

c---
      End If
c---

c-------------------------------
c Inspect the first Inspec columns
c of the linear system
c-------------------------------

      Go to 333   ! comment to display

      Inspec = 6

      write (6,*)
      write (6,*) " prtcl_2d_sys: linear system:"
      write (6,*)

      Do l=1,Ncl2
        write (6,101) (AL(l,m),m=1,Inspec),BL(l)
      End Do

 333  Continue

c----------------------------------

c-----------------------------------------
c Set the size of the Master Linear System
c-----------------------------------------

      Nlin = Ncl2

      write (6,*) " prtcl_2d: sSize of the MLS: ",Nlin

c-----------------------------
c Master linear system defined
c
c Solve by iteration
c-----------------------------

c--------------------------------------------------
c Define MLS particle block sizes (pbs) and indices:
c
c Isz is the size of a particle block
c Imn is the lower index of a particle block
c Imx is the upper index of a particle block
c
c--------------------------------------------------

      Imn(1) = 1              ! first block
      Isz(1) = 2*NE(1)
      Imx(1) = Imn(1)-1+Isz(1)

      Do i=2,Nprtcl            ! subsequent blocks
        Imn(i) = Imx(i-1)+1
        Isz(i) = 2*NE(i)
        Imx(i) = Imn(i)-1+Isz(i)
      End Do

c----------------------------------------
c Define MLS solution block sizes (sbs)
c and indices for iterative solution
c
c Nsz is the size of the solution block
c Nmn is the lower index of a solution block
c Nmx is the upper index of a solution block
c-------------------------------------------

      Ic = 0

      Nmn(1) = 1
      Nsz(1) = 0

      Do j=1,Lump(1)            ! first block
       Ic = Ic+1
       Nsz(1) = Nsz(1)+Isz(Ic)
      End Do
      Nmx(1) = Nmn(1)-1+Nsz(1)

      Do i=2,Nblocks            ! further blocks
       Nmn(i) = Nmx(i-1)+1
       Nsz(i) = 0
       Do j=1,Lump(i)
        Ic = Ic+1
        Nsz(i) = Nsz(i)+Isz(Ic)
       End Do
       Nmx(i) = Nmn(i)-1+Nsz(i)
      End Do

c-------------------------------------------
c Decompose MLS into Nblocks diagonal blocks
c and iterate
c-------------------------------------------

      Iter = 0

   96 Continue       

      Iter = Iter+1

      Do 39 i=1,Nblocks        ! loop over diagonal solution blocks

c------------------------------
c Update the ith block 
c on the right-hand side of MLS
c
c This is done only when the
c decomposition is non-trivial
c------------------------------

        kc = 0 

        Do k=Nmn(i),Nmx(i)     ! run over block rows

         kc = kc+1
         Bblock(kc) = BL(k)
                               ! will skip the diagonal block
                               ! spanning: Nmn-1 < m < Nmx+1
         If(Nmn(i).gt.1) then
           Do m=1,Nmn(i)-1
            Bblock(kc) = Bblock(kc)-AL(k,m)*SOL(m)
           End Do
         End If

         If(Nmx(i).le.Nlin-1) then
          Do m=Nmx(i)+1,Nlin
           Bblock(kc) = Bblock(kc)-AL(k,m)*SOL(m)
          End Do
         End If

        End Do

c-----------------------------------------
c Compute the diagonal block inverses
c and put them in the master inverse matrix
c called: blockinv
c
c Do this only at first iteration
c-----------------------------------------

       If(Iter.gt.1) Go to 73

         kc = 0               ! pick up ith diagonal block 

         Do k=Nmn(i),Nmx(i)

          kc = kc+1
          lc = 0

          Do l=Nmn(i),Nmx(i)
           lc = lc+1
           Ablock(kc,lc) = AL(k,l)
          End Do

         End Do

         Isym_g = 0   ! system is not symmetric
         Iwlpvt = 0   ! pivoting disabled

         call gel_inv    ! compute the inverse
     +
     +      (Nsz(i)
     +      ,Ablock
     +      ,Ainvblock
     +      ,Isym_g
     +      ,Iwlpvt
     +      ,Deter
     +      ,Istop
     +      )

         Do l=1,Nsz(i)
          Do k=1,Nsz(i)
           blockinv(i,l,k) = Ainvblock(l,k)
          End Do
         End Do

c      Do l=1,Nsz(i)
c        write (6,101) (Ablock(l,m),m=1,Nsz(i))
c      End Do

  73  Continue

c------------------------
c      call gauss 
c    +
c    +     (Nsz(i)
c    +     ,Ablock,Bblock,SOLblock
c    +     ,Isym_g
c    +     ,Iwlpvt
c    +     ,Deter
c    +     ,Istop
c    +     )
c--------------------------

c-----------------------------
c solve the ith block equation
c in terms of the inverse
c-----------------------------

       Do l=1,Nsz(i)
        SOLblock(l) = 0.0D0
        Do k=1,Nsz(i)
         SOLblock(l) = SOLblock(l)+blockinv(i,l,k)*Bblock(k)
        End Do
       End Do

c      Do l=1,Nsz(i)
c        write (6,101) (Ablock(l,m),m=1,Nsz(i)),Bblock(l),SOLblock(l)
c      End Do

c----------------------------------
c update the entries of the master
c solution vector 
c corresponding to the ith block
c----------------------------------

       Do l=1,Nsz(i)
        SOL(Nmn(i)+l-1) = SOLblock(l)
       End Do

c--------------
   39  Continue ! end of loop over diagonal solution blocks
c--------------


c-----------------------------------
c save the solution for a posteriori
c error estimates
c-----------------------------------

      Do i=1,Nlin
       SOLsave(i,Iter) = SOL(i)
      End Do

c------------------------
c compute the pointwise error
c and stop iterations if it is
c smaller than tol
c------------------------

      If(Iter.gt.1) then

       Errmax = 0.0

       Do i=1,Nlin
        Diff = abs(SOLsave(i,Iter)-SOLsave(i,Iter-1))
        If(Diff.gt.Errmax) Errmax = Diff
       End Do

       perr(Iter) = Errmax

       write (6,888) Iter,Errmax

       If(Errmax.lt.tol) Go to 98     ! escape

      End If

c---------------------
c Aitken extrapolation
c and pointwise error
c---------------------
c
c     If(Iter.gt.2) then
c       Errmax = 0.
c       Do i=1,Nlin
c        sol0 = SOLsave(i,Iter-2)
c        sol1 = SOLsave(i,Iter-1)
c        sol2 = SOLsave(i,Iter)
c        SOLext(i,Iter) = (sol0*sol2-sol1**2)/(sol2-2.0*sol1+sol0)
c        Diff = abs(SOLext(i,Iter)-SOLext(i,Iter-1))
c        If(Diff.gt.Errmax) Errmax = Diff
c       End Do
c       write (6,*) "Iteration:",Iter," Pointwise extr error : ",Errmax
c       If(Errmax.lt.tol) then
c         Do i=1,Nlin
c          SOL(i) = SOLext(i,Iter)
c         End Do
c         Go to 98
c       End If
c      End If
c
c-------------------------
c Steffensen extrapolation
c-------------------------
c
c     If(Isteff.eq.1.and.Iter.gt.1.and.mod(Iter,2).eq.1) then
c      If(mod(Iter,2).eq.1) then
c       Do i=1,Nlin
c        Sol(i) = SOLext(i,Iter) 
c       End Do
c      End If
c     End If
c----------------------------------------------------

c---
      If(Nblocks.eq.1) Go to 98   ! MLS solved at once
c---
      If(Iter.lt.Niter) Go to 96
c---

      write (6,*) " prtcl_2d_sys: Maximum number of iterations "
      write (6,*) "               exceeded; will abort"
      write (6,*)

      Istop = 1

c--------------------

  98  Continue

c     write (6,*)
c     write (6,*) " Solution of the MLS"
c     write (6,*)
c     write (6,111) (SOL(l),l=1,Nlin)

c----------------------------
c compute rate of convergence
c a posteriori
c----------------------------

      If(Iter.gt.1) then

      open (15,file="prtcl_2d.error")
 
      write (6,*) 
      write (6,*) " Error in iterations"
      write (6,*)

      write (15,*) Iter-1," Global and pointwise error in iterations"

      Do i=1,Iter-1

       error(i) = 0.0D0

       Do j=1,Nlin
         error(i) = error(i)+(SOLsave(j,i)-SOL(j))**2
       End Do

       error(i) = Dsqrt(error(i))/Nlin
       error(i) = log10(error(i))

       write (06,103) i,error(i),perr(i)
       write (15,103) i,error(i),perr(i)

      End Do

c------------------------------------------------
c     write (6,*) 
c     write (6,*) " Extrapolated error iterations"
c     write (6,*)
c
c     write (15,*) Iter-2," Extrapolated error"
c
c     Do i=2,Iter-1
c      error(i) = 0.0
c      Do j=1,Nlin
c        error(i) = error(i)+(SOLext(j,i)-SOL(j))**2
c      End Do
c      error(i) = sqrt(error(i))/Nlin
c      error(i) = log10(error(i))
c      write (06,103) i,error(i)
c      write (15,103) i,error(i)
c      End Do
c------------------------------------------------

      write (15,103) Null
      close (15)

      End If

c---------------------------------
c  Done with the solution
c
c  Distribute to traction matrices
c---------------------------------

      Ib = 0        ! block counter

      Do i=1,Nprtcl

       Do j=1,NE(i)
        Ibj  = Ib+j
        Ibji = Ibj + NE(i)
        fx(i,j) = visc*SOL(Ibj)
        fy(i,j) = visc*SOL(Ibji)
       End Do

       Ib = Ib+2*NE(i)

c      If(Ireg.eq.1) then     ! system was reduced
c        Ib = Ib-1
c        fy(I,NE(i)) = 0.0
c      End If

      End Do

c-----
c Done
c-----

 101  Format (80(1x,f7.3))
 103  Format (1x,i3,5(1x,f15.10))
 104  Format (1x,i3,20(1x,f9.5))
 111  Format (80(10(1x,f7.3),/))
 112  Format (80(1x,f10.7))

 888  Format (1x,"Iter: ",i3," Point corr: ",f15.10)

      Return
      End
