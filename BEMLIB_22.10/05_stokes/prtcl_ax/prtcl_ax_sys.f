      subroutine prtcl_ax_sys
     +
     +  (Iflow
     +  ,Iprec
     +  ,Ireduce
     +  ,Nblocks
     +  ,Lump
     +  ,Niter
     +  ,tol
     +  ,Istop
     +  )

c-----------------------------------------
c FDLIB and BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c---------------------------------
c Prepare and solve the 
c Master Linear system (MLS)
c---------------------------------

c------------------------------------------
c Capacity:
c
c   25 particles
c   64 elements along each particle
c
c LEGEND:
c -------
c
c Imn(i): Low  index of ith particle in master matrix
c Imx(i): High index of ith particle in master matrix
c Isz(i):       size of ith particle in master matrix
c
c Ncl:    number of collocation points
c
c------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    NE(25),  Itp(25)
      Dimension xcntr(25),ycntr(25)
      Dimension axis1(25),axis2(25),tilt(25)

      Dimension  xw(25,65),yw(25,65),tw(25,65)
      Dimension  fx(25,64),fy(25,64)

      Dimension   X0(3200),  Y0(3200),T0(3200)
      Dimension  ux0(3200), uy0(3200)
      Dimension elar(3200)
      Dimension vnX0(3200),vnY0(3200)
      Dimension tnX0(3200),tnY0(3200)

c---
c related to the linear solver
c---

      Dimension Imn(25),Imx(25),Isz(25)
      Dimension Nmn(25),Nmx(25),Nsz(25),Lump(25)

      Dimension AL(3200,3200),BL(3200),SOL(3200) 

      Dimension Ablock(500,500),Ainvblock(500,500)
      Dimension Bblock(500),SOLblock(500) 

      Dimension blockinv(25,128,128) 

c---
c related to iterations
c---

      Dimension SOLsave(1600,200)

      Dimension error(200)    ! rms error
      Dimension perr(200)     ! pointwise error

c     Dimension SOLext(1600,200)

c--------------
c common blocks
c--------------

      common/INTGR1/Nprtcl,NGL,NE,Itp,Ncl

      common/points/xw,yw,tw
      common/particles/xcntr,ycntr,axis1,axis2,tilt

      common/colloc1/fx,fy,ux0,uy0
      common/colloc2/vnx0,vny0,elar
      common/colloc3/tnx0,tny0
      common/colloc4/x0,y0,t0

      common/REAL1/visc,Uinf,wall,pg,RL,sc

      common/matrix/AL,BL

      common/piii/pi,pi2,pi4,pi8

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
c------------------------------

      Ic = 0       ! collocation point counter
      Ib = 0       ! block counter

      Do 31 i=1,Nprtcl     ! loop over particles

       Do 32 j=1,NE(i)     ! loop over elements

       Ic = Ic+1

c      write (6,*) " prtcl_ax_sys: collocating at point :",Ic

       Ibj  = Ib  + j      ! row entry for x component of BIE
       Ibji = Ibj + NE(i)  ! row entry for y component of BIE

c---
c right-hand side
c---

       BL(Ibj ) = - pi8*ux0(Ic) 
       BL(Ibji) = - pi8*uy0(Ic)

c      sum1 = 0.0   ! for debugging
c      sum2 = 0.0
c      sum3 = 0.0
c      sum4 = 0.0
c      sum5 = 0.0
c      sum6 = 0.0

       Kc = 0        ! element counter
       kb = 0        ! block counter

       Do 33 k=1,Nprtcl       ! loop over particles
        Do 34 l=1,NE(k)       ! loop over elements

         Kc = Kc+1

         kbl  = kb  + l         ! column entry for x traction
         kblk = kbl + ne(k)     ! column entry for y traction

         X1 = XW(K,L)          ! first element point
         Y1 = YW(K,L)
         T1 = TW(K,L)

         X2 = XW(K,L+1)        ! second element point
         Y2 = YW(K,L+1)
         T2 = TW(K,L+1)

         Ising = 0
         If(Ic.eq.Kc) Ising = 1

         call prtcl_ax_slp
     +
     +      (Iflow
     +      ,X0(Ic),Y0(Ic),t0(Ic)
     +      ,X1,Y1,T1
     +      ,X2,Y2,T2
     +      ,NGL
     +      ,Ising
     +      ,Itp(k)
     +      ,xcntr(k),ycntr(k)
     +      ,axis1(k),axis2(k),tilt(k)
     +      ,QQxx,QQxy
     +      ,QQyx,QQyy
     +      )

         AL(Ibj, kbl)  = QQxx
         AL(Ibj, kblk) = QQxy
         AL(Ibji,kbl)  = QQyx
         AL(Ibji,kblk) = QQyy

c        write (6,*) QQxx
c        pause

c        sum5 = sum5 + AL(Ibj, Kbl) *vnx0(Kc)   ! for debugging
c    +               + AL(Ibj, kblk)*vny0(Kc)
c        sum6 = sum6 + AL(Ibji,Kbl) *vnx0(Kc)
c    +               + AL(Ibji,kblk)*vny0(Kc)

   34    Continue

         Kb = Kb + 2*NE(K)     ! shift by one particle block

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

      Ib = Ib + 2*NE(i)     ! shift by one particle block

  31  Continue

c-------------------------------
c Inspect the first inspec columns
c of the linear system
c-------------------------------
c
c     inspec = 5
c     write (6,*)
c     Do l=1,Ncl2
c       write (6,101) (AL(l,m),m=1,inspec),BL(l)
c     End Do
c     pause
c----------------------------------

c-------------------------
c precondition and reduce
c-------------------------

      If(Iprec.eq.1)   call precondition ()
      If(Ireduce.eq.1) call reduce ()

c-------------------------------
c Inspect the first inspec columns
c of the linear system
c-------------------------------
c
c     inspec = 5
c     write (6,*)
c     Do l=1,Ncl2
c       write (6,101) (AL(l,m),m=1,inspec),BL(l)
c     End Do
c----------------------------------

c-----------------------------------------
c Set the size of the Master Linear System
c
c If reduction was done, system size is reduced 
c by the number of particles
c-----------------------------------------

      Nlin = Ncl2

      write (6,*) " prtcl_ax_sys: size of the linear system:",Nlin

c--------------------------------------------------
c Define MLS particle block sizes (pbs) and indices
c
c Isz is the size of a particle block
c Imn is the lower index of a particle block
c Imx is the upper index of a particle block
c--------------------------------------------------

      Imn(1) = 1              ! first particle block
      Isz(1) = 2*NE(1)
      Imx(1) = Imn(1)+Isz(1)-1

      Do i=2,Nprtcl            ! subsequent particle blocks
        Imn(i) = Imx(i-1)+1
        Isz(i) = 2*NE(i)
        Imx(i) = Imn(i)+Isz(i)-1
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

      Do j=1,Lump(1)
       Ic = Ic+1
       Nsz(1) = Nsz(1)+Isz(Ic)
      End Do
      Nmx(1) = Nmn(1)+Nsz(1)-1

      Do i=2,Nblocks
       Nmn(i) = Nmx(i-1)+1
       Nsz(i) = 0
       Do j=1,Lump(i)
        Ic = Ic+1
        Nsz(i) = Nsz(i)+Isz(Ic)
       End Do
       Nmx(i) = Nmn(i)+Nsz(i)-1
      End Do

c-------------------------------------------
c Decompose MLS into Nblocks diagonal blocks
c and iterate
c-------------------------------------------

      Iter = 0

   96 Continue       

      Iter = Iter+1

      Do i=1,Nblocks        ! loop over diagonal solution blocks

c------------------------------
c Update ith block 
c on the right-hand side of MLS
c
c This is done only when the
c decomposition is non-trivial
c------------------------------

        kc = 0 

        Do k=Nmn(i),Nmx(i)

         kc = kc+1
         Bblock(kc) = BL(k)

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

         call gauss_inv                   ! compute the inverse
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
c    +            (Nsz(i)
c    +            ,Ablock,Bblock,SOLblock
c    +            ,Isym_g
c    +            ,Iwlpvt
c    +            ,Deter
c    +            ,Istop
c    +            )
c--------------------------

c-----------------------------
c solve the ith block equation
c in terms of the inverse
c-----------------------------

       Do l=1,Nsz(i)
        SOLblock(l) = 0.0
        Do k=1,Nsz(i)
         SOLblock(l) = SOLblock(l)+blockinv(i,l,k)*Bblock(k)
        End Do
       End Do

c      Do l=1,Nsz(i)
c        write (6,101) (Ablock(l,m),m=1,Nsz(i)),Bblock(l),SOLblock(l)
c      End Do

c----------------------------------
c update the partition of the master
c solution  vector 
c corresponding to the ith block
c----------------------------------

       Do l=1,Nsz(i)
        SOL(Nmn(i)+l-1) = SOLblock(l)
       End Do

      End Do       ! end of loop over diagonal solution blocks

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

      write (6,*) " prtcl_ax_sys: maximum number of iterations exceeded"
      write (6,*) "               aborting"

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
c
c     If(Iter.gt.1) then
c
c     open (15,file="prtcl_ax.error")
c 
c     write (6,*) 
c     write (6,*) " prtcl_ax_sys: Error in iterations"
c     write (6,*)
c 
c     write (15,*) Iter-1," Global and pointwise error in iterations"
c 
c     Do i=1,Iter-1
c
c      error(i) = 0.0D0
c
c      Do j=1,Nlin
c        error(i) = error(i)+(SOLsave(j,i)-SOL(j))**2
c      End Do
c
c      error(i) = sqrt(error(i))/Nlin
c      error(i) = log10(error(i))
c
c      write (06,103) i,error(i),perr(i)
c      write (15,103) i,error(i),perr(i)
c
c     End Do

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
c     write (15,103) Null
c
c     close (15)
c
c     End If
c----------------

c----------------------------------------------
c  Distribute the solution to traction matrices
c-----------------------------------------------

      Ib = 0        ! block counter

      Do i=1,Nprtcl

       Do j=1,NE(i)
        Ibj  = Ib+j
        Ibji = Ibj + NE(i)
        fx(i,j) = visc*SOL(Ibj)
        fy(i,j) = visc*SOL(Ibji)
       End Do

       Ib = Ib+2*NE(i)

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
