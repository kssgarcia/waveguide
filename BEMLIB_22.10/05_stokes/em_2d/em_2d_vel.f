      subroutine em_2d_vel
     +
     +  (IS_slp
     +  ,IS_dlp
     +  ,Isolve,JGS
     +  ,eps,Nter
     +  ,Iflow
     +  ,Idfl
     +  ,Istop
     +  )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c----------------------------------------
c Solve an integral equation for the velocity
c over the interfaces
c
c CAPACITY:
c
c When all viscosities are equal:
c
c 49 interfaces
c 64 points around each interface (66 including wrap around)
c
c When viscosities are different:
c
c 05 interfaces
c 64 points around each interface (66 including wrap around)
c
c
c Global variables:
c
c the g in the end indicates "global"
c-----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension xg(0:066,49),yg(0:066,49)
      Dimension ug(0:066,49),vg(0:066,49)

      Dimension  XCg(066,49), YCg(066,49),  Rg(066,49),   sg(066,49)
      Dimension TH1g(066,49),TH2g(066,49),TH3g(066,49),orntg(066,49)

      Dimension Uinfg(066,49),Vinfg(066,49)
      Dimension Uslpg(066,49),Vslpg(066,49)
      Dimension Unewg(066,49),Vnewg(066,49)

      Dimension tpxx(066,49),tpyx(066,49)
      Dimension tpxy(066,49),tpyy(066,49)

      Dimension DLMxx(066,05,066,05),DLMxy(066,05,066,05)
      Dimension DLMyx(066,05,066,05),DLMyy(066,05,066,05)

c---
c properties of the drops, interfaces
c---

      Dimension  NSG(49),NSG1(49),NSG2(49)
      Dimension  vsd(49),Ivsg(49)
      Dimension rhod(49)
      Dimension elst(49)

      Dimension CF1(49),CF2(49),CF3(49),CF4(49)
      Dimension xcnt(49),ycnt(49)

c---
c local variables
c---

      Dimension X(0:200),Y(0:200)
      Dimension U(0:200),V(0:200)

      Dimension  XC(200), YC(200),  R(200),   S(200)
      Dimension TH1(200),TH2(200),TH3(200),ORNT(200)

c---
c Matrices and vectors for direct solution
c system size is Pro ( 2*NSG(i))
c naximum size is 25 * 2 * 64 = 3200
c---

      Dimension amat(3200,3200),rhs(3200),sln(3200)

c-----------------------------------
c common blocks for global variables
c-----------------------------------

      common/XXYYg/Xg,Yg
      common/ARCCg/XCg,YCg,Rg,Sg,TH1g,TH2g,TH3g,ORNTg
      common/UUVVg/Ug,Vg

c---
c common blocks for local variables
c---

      common/XXYY/X,Y
      common/ARCC/XC,YC,R,S,TH1,TH2,TH3,ORNT
      common/UUVV/U,V

c---
c various common blocks
c---

      common/ancR1/gx,gy
      common/ancR2/vs1,vsd
      common/ancR4/rho1,rhod,elst
      common/ancR6/shrt,U1,U2,pg

      common/hhhh/h,hh,h2,h3,h4,hs
      common/teihos/wall

      common/ancI1/Ndrops,NSG,NSG1,NSG2
      common/ancI4/Ivs,Ivsg

      common/piii/pi,pih,piq,pi2,pi4,pi6,pi8,srpi

c--------------
c incident flow
c--------------

      Do j=1,Ndrops
       Do i=1,NSG1(j)

c------------------------
       if(Iflow.eq.1) then     ! straining flow
c------------------------

         Uinfg(i,j) =  shrt*xg(i,j)
         Vinfg(i,j) = -shrt*yg(i,j)

c---------------------------------------
       else if(Iflow.eq.2.or.Iflow.eq.5) then  !  simple shear flow
c---------------------------------------

         Uinfg(i,j) = shrt*yg(i,j)
         Vinfg(i,j) = 0.0D0

c----------------------------------------------
       else if(Iflow.eq.10.or.Iflow.eq.11) then  ! shear flow above a wall
c----------------------------------------------

         Uinfg(i,j) = shrt*(yg(i,j)-wall)
         Vinfg(i,j) = 0.0D0

c------------------------------
      else if(Iflow.eq.20) then     ! periodic channel flow
c------------------------------

         Uinfg(i,j) = ( U1*(h-yg(i,j))+U2*(h+yg(i,j)) )/(2.0D0*h)
     +               + 0.5D0*pg*(hs-yg(i,j)**2)
     +               + (0.5D0*gx*rho1/vs1) *yg(i,j)*(h2-yg(i,j))
         Vinfg(i,j) = 0.0D0

c------------
      else
c------------

         write (6,*) " em_2d_vel: this flow cannot be selected"
         stop

c------------
      end if
c------------

       End Do
      End Do

c----------------------------------
c compute the single-layer potential
c----------------------------------

      write (6,*) " em_2d_vel: computing the SLP"

      call em_2d_slp 
     +
     +   (IS_slp
     +   ,Iflow
     +   ,uslpg,vslpg
     +   ,Istop
     +   )

      write (6,*) " em_2d_vel: SLP computed"

c---
c inspect slp
c---

c      Do j=1,Ndrops
c       write (6,*)
c       Do i=1,NSG(j)
c        write (6,100) j,i,uslpg(i,j),vslpg(i,j)
c       End Do
c      End Do

c--------------------------------------------------

c---
c all viscosity ratios are equal to unity
c---

      if(Ivs.eq.0) then

       Do j=1,Ndrops
        Do i=1,NSG(j)
         ug(i,j) = Uinfg(i,j)+uslpg(i,j)/vs1
         vg(i,j) = Vinfg(i,j)+vslpg(i,j)/vs1
c        write (6,*) j,i,Uinfg(i,j),Vinfg(i,j),uslpg(i,j),vslpg(i,j)
        End Do
       End Do

      Go to 99

      end if

c-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

c---
c viscosity ratios different than unity
c---

c-------------------------
c Generate the dlp matrix
c by the method of impulses
c-------------------------

      write (6,*) " em_2d_vel: generating the dlp-matrix"

      Do 1 l=1,Ndrops       ! run over drops

c---
      if(Ivsg(l).eq.0) then        ! of lamda(l) = 1 bypass

        Do j=1,NSG(l)
         Do k=1,Ndrops
          Do i=1,NSG(k)
           DLMxx(i,k,j,l) = 0.0D0
           DLMyx(i,k,j,l) = 0.0D0
           DLMxy(i,k,j,l) = 0.0D0
           DLMyy(i,k,j,l) = 0.0D0
          End Do
         End Do
        End Do
        Go to 1

      end if
c---

      Do j=1,NSG(l)

          XCS =   XCg(j,l)
          YCS =   YCg(j,l)
         RADS =    Rg(j,l)
         TH1S =  TH1g(j,l)
         TH2S =  TH2g(j,l)
         TH3S =  TH3g(j,l)
        ORNTS = ORNTg(j,l)

        call em_2d_dlp 
     +
     +    (l,j
     +    ,XCS,YCS,RADS
     +    ,TH1S,TH2S,TH3S,ORNTS
     +    ,IS_dlp
     +    ,Iflow
     +    ,Tpxx,Tpyx
     +    ,Tpxy,Tpyy
     +    ,Istop
     +    )

         Do k=1,Ndrops
          Do i=1,NSG(k)
            DLMxx(i,k,j,l) = tpxx(i,k)
            DLMyx(i,k,j,l) = tpyx(i,k)
            DLMxy(i,k,j,l) = tpxy(i,k)
            DLMyy(i,k,j,l) = tpyy(i,k)
          End Do
         End Do

      End Do

c---

  1   Continue

c---------------------------------------
c Combine the incident flow with the slp
c---------------------------------------

      Do m=1,Ndrops
        cf = 1.0D0/(1.0D0+vsd(m)/vs1)
        Do i=1,NSG(m)
c        write (6,100) m,i,Uslpg(i,m),Uinfg(i,m),vslpg(i,m),Vinfg(i,m)
         uslpg(i,m) = cf*(uinfg(i,m)+uslpg(i,m)/vs1)
         vslpg(i,m) = cf*(vinfg(i,m)+vslpg(i,m)/vs1)
        End Do
      End Do

c-----------------------------
c Prepare the iteration matrix
c-----------------------------

      Do m=1,Ndrops
       Do l=1,Ndrops
        cf = (vs1-vsd(l))/(vs1+vsd(m))
        cf = cf/pi2
        Do i=1,NSG(m)
          Do j=1,NSG(l)
            DLMxx(i,m,j,l) = cf*DLMxx(i,m,j,l)
            DLMyx(i,m,j,l) = cf*DLMyx(i,m,j,l)
            DLMxy(i,m,j,l) = cf*DLMxy(i,m,j,l)
            DLMyy(i,m,j,l) = cf*DLMyy(i,m,j,l)
          End Do
        End Do
       End Do
      End Do

c----------------
c DIRECT SOLUTION
c----------------

      if(Isolve.eq.1) then

      write (6,*) " Solving the linear system "

      mats = 0          ! matrix size
      Do k=1,Ndrops
        mats = mats+NSG(k)
      End Do
      mats = 2*mats

c---
c set up a linear system
c---
 
      Imarkk = 0     ! defines row blocks

      Do k=1,Ndrops
 
        Do i=1,NSG(k)

        ii1 = Imarkk+i
        ii2 = ii1+NSG(k)

          Imarkl = 0            ! defines column blocks

          Do l=1,Ndrops
            Do j=1,NSG(l)
             jj1 = Imarkl+j
             jj2 = jj1+NSG(l)
             amat(ii1,jj1) = - DLMxx(i,k,j,l)
             amat(ii2,jj1) = - DLMxy(i,k,j,l)
             amat(ii1,jj2) = - DLMyx(i,k,j,l)
             amat(ii2,jj2) = - DLMyy(i,k,j,l)
            End Do
            Imarkl = Imarkl+2*NSG(l)
          End Do

         amat(ii1,ii1) = amat(ii1,ii1)+1.0 ! left-hand side
         amat(ii2,ii2) = amat(ii2,ii2)+1.0 ! left-hand side
         rhs (ii1)     = uslpg(i,k)
         rhs (ii2)     = vslpg(i,k)

         End Do

       Imarkk = Imarkk+2*NSG(k)

      End Do
 
c---
c Gauss elimination
c---

      Iwlpvt = 1      ! will pivot
      Isym_g = 0      ! coefficient matrix is not symmetric

      call gauss_large 
     +
     +  (mats,amat,rhs,sln
     +  ,Isym_g,Iwlpvt
     +  ,l,u
     +  ,det
     +  ,Istop
     +  )

c---
c distribute the  solution
c---

       Imarkk = 0

       Do k=1,Ndrops
          Do i=1,NSG(k)
           ii1 = Imarkk+i
           ii2 = ii1+NSG(k)
c          ug(i,k) = sln(ii1)
c          vg(i,k) = sln(ii2)
         End Do
          Imarkk = Imarkk+2*NSG(k)
       End Do
 
c---
c Done with the direct solution
c---

        write (6,*) " caps_2d_vel: linear system  solved"
 
        Go to 99

      end if

c-------------------
c NEUMANN ITERATIONS
c-------------------

      write (6,*) " caps_2d_vel: begin iterations"

      Iter  = 0

  98  Continue

      Iter  = Iter + 1
      Iflag = 0                 ! flag for convergence

c------------------------------------------------------
c  Deflation
c  ---------
c
c  Idfl = 0 no deflation
c  Idfl = 1 one eigenfunction  will be removed 
c  Idfl = 2 all four eigenfunctions will be removed 
c-----------

                             ! four modes of deflation
      Do j=1,Ndrops
       CF1(j) = 0.0          ! normal vector
       CF2(j) = 0.0          ! translation along the x axis
       CF3(j) = 0.0          ! translation along the y axis
       CF4(j) = 0.0          ! rotation about the z axis
      End Do

      if(Idfl.ne.0) then

c---
c  Transfer into local vectors
c---

        Do j=1,Ndrops

          rkap = (vs1-vsd(j))/(vs1+vsd(j))
 
          Do i=1,NSG1(j)
              x(i) =    xg(i,j)
              y(i) =    yg(i,j)
             xc(i) =   xcg(i,j)
             yc(i) =   ycg(i,j)
              r(i) =    rg(i,j)
              s(i) =    sg(i,j)
            th1(i) =  th1g(i,j)
            th2(i) =  th2g(i,j)
            th3(i) =  th3g(i,j)
           ornt(i) = orntg(i,j)
              u(i) =    ug(i,j)
              v(i) =    vg(i,j)
          End Do

          call deflation 
     +
     +      (Idfl
     +      ,NSG(j)
     +      ,U,V
     +      ,CF1(j),CF2(j),CF3(j),CF4(j)
     +      ,xcnt(j),ycnt(j)
     +      )

          CF1(j) = rkap * CF1(j)
          CF2(j) = rkap * CF2(j)
          CF3(j) = rkap * CF3(j)
          CF4(j) = rkap * CF4(j)

        End Do

      end if

c---------------------------------------------------------

c-----------
c  Updating
c-----------

      corr = 0.0

      Do 5 m=1,Ndrops
      Do 5 i=1,NSG(m)

        Udl = 0.0D0
        Vdl = 0.0D0

        Do l=1,Ndrops

c       sum1 = 0.0    ! for debugging
c       sum2 = 0.0
c       sum3 = 0.0
c       sum4 = 0.0

          Do j=1,NSG(l)
            Udl = Udl + Ug(j,l)*DLMxx(i,m,j,l)
     +                + Vg(j,l)*DLMyx(i,m,j,l)
            Vdl = Vdl + Ug(j,l)*DLMxy(i,m,j,l)
     +                + Vg(j,l)*DLMyy(i,m,j,l)
c           sum1 = sum1+DLMxx(i,m,j,l)
c           sum2 = sum2+DLMyx(i,m,j,l)
c           sum3 = sum3+DLMxy(i,m,j,l)
c           sum4 = sum4+DLMyy(i,m,j,l)
          End Do

c       write (6,101) sum1,sum2,sum3,sum4

        End Do

        Unewg(i,m) = uslpg(i,m)+Udl
        Vnewg(i,m) = vslpg(i,m)+Vdl

c---
c deflation
c---
  
        Unewg(i,m) = Unewg(i,m)
     +             - CF1(m) * cos(th2g(i,m))*orntg(i,m)
     +             + CF2(m)
     +             + CF4(m)*(-yg(i,m)+ycnt(m))

        Vnewg(i,m) = Vnewg(i,m)
     +             - CF1(m) * sin(th2g(i,m))*orntg(i,m)
     +             + CF3(m)
     +             + CF4(m)*(xg(i,m)-xcnt(m))

c---
c increments
c---

        DU = Unewg(i,m)-Ug(i,m)
        DV = Vnewg(i,m)-Vg(i,m)

        if(JGS.eq.2) then         ! Gauss-Siedel
          Ug(i,m) = Unewg(i,m)
          Vg(i,m) = Vnewg(i,m)
        end if

        if(Dabs(DU).ge.eps.or.Dabs(DV).ge.eps) Iflag = 1 
        if(Dabs(DU).ge.Corr) Corr = abs(DU)
        if(Dabs(DV).ge.Corr) Corr = abs(DV)

c       write (6,199) i,DU,DV,UDL,VDL

  5   Continue

c---
c End of updating
c---

c------------
c diagnostics
c------------

      Do j=1,Ndrops
       if(Idfl.eq.0) write (6,195) Iter,Corr
       if(Idfl.eq.1) write (6,196) Iter,Corr,CF1(j)
       if(Idfl.eq.2) write (6,197) Iter,Corr,CF1(j),CF2(j)
     +                                      ,CF3(j),CF4(j)
      End Do

c-------------------------------------

      if(JGS.eq.1) then         ! Jacobi

        Do m=1,Ndrops
          Do i = 1,NSG(m)
            Ug(i,m) = Unewg(i,m)
            Vg(i,m) = Vnewg(i,m)
          End Do
          Ug(nsg1(m),m) = Ug(1,m)   ! needed for deflation
          Vg(nsg1(m),m) = Vg(1,m)   ! needed for deflation
        End Do

      end if

c--------------------------------------

      if(Iter.gt.Nter) then
        write (6,*) " Iterations failed to converge"
        Go to 91
      end if

      if(Iflag.eq.1) Go to 98 

c-------------------------
c Done with the iterations
c-------------------------

c--------------------------------

  91  Continue

      if(Isolve.eq.1) Idfl = 0

c-------------------
c  if deflation of RBM was done
c  rectify the solution
c-------------------

      if(Idfl.eq.2) then

       Do j=1,Ndrops
        rkap = (vs1-vsd(j))/(vs1+vsd(j))
        cf = 1.0/(1.0+rkap)
        Do i =1,NSG(j)
         Ug(i,j) = ug(i,j)-cf*(CF2(j)-CF4(j)*(yg(i,j)-ycnt(j)))
         Vg(i,j) = vg(i,j)-cf*(CF3(j)+CF4(j)*(xg(i,j)-xcnt(j)))
        End Do
       End Do

      end if

c-----------------------------------------------

  99  Continue

c------------
c wrap around
c------------

      Do j=1,Ndrops

        Ug(0,j) = Ug(NSG(j),j)
        Vg(0,j) = Vg(NSG(j),j)

        Ug(NSG1(j),j) = Ug(1,j)
        Vg(NSG1(j),j) = Vg(1,j)

        Ug(NSG2(j),j) = Ug(2,j)
        Vg(NSG2(j),j) = Vg(2,j)

      End Do

c-----
c done
c-----

 100  Format (1X,I2,1X,I2,9(1X,F10.5))
 101  Format (9(1X,F10.5))

 195  Format (' Iter:',I3,' Max Corr=',F9.6)
 196  Format (' Iter:',I3,' Max Corr=',F9.6,' CF1= ',F7.5)
 197  Format (' Iter:',I3,' Max Corr=',F9.6,' CF1= ',F7.5,
     +          ' CF2 = ',f7.5," CF3 =",F7.5," CF4=",F7.5)
 198  Format (' Evaluating the SLP at point: ',I2)
 199  Format (' Increments at point ',I2,4(F10.5))

      return
      end
