      subroutine caps_2d_vel
     +
     +  (IS_slp
     +  ,IS_dlp
     +  ,Isolve
     +  ,JGS
     +  ,Iflow
     +  ,Isym
     +  ,Idfl
     +  ,roexp,mexp,mexp_foam
     +  ,Istop
     +  )

c=======================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=======================================

c---------------------------------------------
c Solve an integral equation over an interface
c and return the velocity
c---------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension X(0:900),Y(0:900)
      Dimension U(0:900),V(0:900)

      Dimension XC (900),YC (900),  R(900),   S(900)
      Dimension TH1(900),TH2(900),TH3(900),ORNT(900)

      Dimension Uinf(900),Vinf(900)
      Dimension Uslp(900),Vslp(900)
      Dimension Unew(900),Vnew(900)

      Dimension tpxx(900),tpyx(900)
      Dimension tpxy(900),tpyy(900)

      Dimension DLMxx(500,500),DLMxy(500,500)
      Dimension DLMyx(500,500),DLMyy(500,500)

      Dimension Amat(500,500),rhs(500),sln(500)

c---
c eigenvector of the influence matrix:
c---

      Dimension evt(500)

c--------------
c common blocks
c--------------

      common/XXYY/X,Y
      common/ARCC/XC,YC,R,S,TH1,TH2,TH3,ORNT
      common/UUVV/U,V

      common/ancR1/thet0,gx,gy
      common/ancR2/vs1,vs2,vsr,vsr1,vsr2,rkap,rkap2,eps
      common/ancR4/Drho
      common/ancR6/shrt,yref,pg,ant_c1,ant_c2

      common/hhhh/h,hh,h2,h3,h4,hs
      common/teihos/wall

      common/ancI1/NSG,NSG1,NSG2,NSGM,NSGQ
      common/ancI4/Ivs,Nter

      common/piii/pi,pih,piq,pi2,pi4,pi6,pi8,srpi

      common/aaaa/a11,a12,a21,a22
      common/ewew/ew,tau

c---------
c announce
c---------

      write (6,*) " caps_2d_vel: Solving the integral equation"

c-----------------------------
c  Compute velocity at points:
c  1, ..., Imax
c-----------------------------

      if(Isym.eq.0) then

        Imax = NSG

      else if(Isym.eq.2) then

        NSGQ = NSG/4+1
        NSGM = NSG/2+1
        NSGT = NSGM+NSGQ-1
        Imax = NSGQ

      end  if

      Imax2 = 2*Imax

c--------------
c incident flow
c--------------

      Do i=1,Imax
 
c-------------------------
       if(Iflow.eq.1) then     ! straining flow
c-------------------------

         Xs = X(i)**2
         Ys = Y(i)**2
         antx = 1.0D0+2.0D0*ant_c1*(Xs-3.0D0*Ys)+
     +                      ant_c2*(Xs+3.0D0*Ys)
         anty = 1.0D0+2.0D0*ant_c1*(3.0D0*Xs-Ys)+
     +                      ant_c2*(3.0D0*Xs+Ys)
         Uinf(i) =  X(i)*antx*shrt
         Vinf(i) = -Y(i)*anty*shrt

c--------------------------------------------
       else if(Iflow.eq.2.or.Iflow.eq.5) then  ! simple shear flow
c--------------------------------------------

         Uinf(i) = shrt*Y(i)
         Vinf(i) = 0.0D0

c-------------------------
       else if(Iflow.eq.6) then   ! doubly-periodic expanding flow
c-------------------------

         Uinf(i) = 0.0D0
         Vinf(i) = 0.0D0

c-------------------------
       else if(Iflow.eq.10.or.Iflow.eq.11) then  ! shear flow above a wall
c-------------------------

         Uinf(i) = shrt*(Y(i)-wall)
         Vinf(i) = 0.0D0

c-------------------------
       else if(Iflow.eq.20) then    ! periodic collection in channel Couette
c-------------------------

         Uinf(i) = shrt*(Y(i)-yref)
     +            + 0.5D0*pg/vs1*(hs-y(i)**2)
     +            + (0.5D0*gx*rho1/vs1) * y(i)*(hs-y(i)**2) 
         Vinf(i) = 0.0D0

c----------
       else
c----------

         write (6,*) " caps_2d_vel: this flow cannot be selected"
         stop

c------------
       end if
c------------

      End Do

c-----------------------------------
c Compute the single-layer potential
c-----------------------------------

      write (6,*) " caps_2d_vel: computing the SLP"

      call caps_2d_slp 
     +
     +   (IS_slp
     +   ,Isym,Imax
     +   ,Iflow
     +   ,Uslp,Vslp
     +   ,Istop
     +   )

      write (6,*) " caps_2d_vel: SLP computed" 

c---
c look at the slp
c---

c      Do i = 1,Imax
c       write (6,100) j,i,Uslp(i),Vslp(i)
c      End Do

c--------------------------------------------------

c---
c viscosity ratio equal to 1
c
c Add the incident flow and return
c---

      if(Ivs.eq.0) then

       Do i=1,Imax
         U(i) = Uinf(i)+Uslp(i)/vs1
         V(i) = Vinf(i)+Vslp(i)/vs1
       End Do

      Go to 99

      end if

c=======================================

c---
c viscosity ratio different than unity
c---

c-------------------------
c Generate the DLP matrix
c by the method of impulses
c-------------------------

      write (6,*) " caps_2d_vel: Generating the dlp-matrix"

      Do j=1,NSG

        XCS   =   XC(j)
        YCS   =   YC(j)
        RADS  =    R(j)
        TH1S  =  TH1(j)
        TH2S  =  TH2(j)
        TH3S  =  TH3(j)
        ORNTS = ORNT(j)

        call caps_2d_dlp 
     +
     +   (j
     +   ,XCS,YCS,RADS
     +   ,TH1S,TH2S,TH3S,ORNTS
     +   ,IS_dlp
     +   ,Imax
     +   ,Iflow
     +   ,Tpxx,Tpyx,Tpxy,Tpyy
     +   ,Istop
     +   )

        Do i=1,Imax
          DLMxx(i,j) = tpxx(i)
          DLMyx(i,j) = tpyx(i)
          DLMxy(i,j) = tpxy(i)
          DLMyy(i,j) = tpyy(i)
c         write (6,100) i,DLMxx(i,j),DLMyx(i,j)
c    +                   ,DLMxy(i,j),DLMyy(i,j)
        End Do

      End Do

c--------------------
c If interface is symmetric with respect to
c x=0 and y = 0
c contract the dlp-matrix
c--------------------

      if(Isym.eq.2) then

        Do i=1,Imax
         DLMxx(i,1) = DLMxx(i,1)-DLMxx(i,NSGM)
         DLMyx(i,1) = DLMyx(i,1)-DLMyx(i,NSGM)
         DLMxy(i,1) = DLMxy(i,1)-DLMxy(i,NSGM)
         DLMyy(i,1) = DLMyy(i,1)-DLMyy(i,NSGM)
        End Do

        Do j=2,NSGQ-1
         j1 = NSGM-j+1
         j2 = NSGM+j-1
         j3 = NSG -j+2
         Do i=1,Imax
          DLMxx(i,j) = DLMxx(i,j)-DLMxx(i,j1)-DLMxx(i,j2)+DLMxx(i,j3)
          DLMxy(i,j) = DLMxy(i,j)-DLMxy(i,j1)-DLMxy(i,j2)+DLMxy(i,j3)
          DLMyx(i,j) = DLMyx(i,j)+DLMyx(i,j1)-DLMyx(i,j2)-DLMyx(i,j3)
          DLMyy(i,j) = DLMyy(i,j)+DLMyy(i,j1)-DLMyy(i,j2)-DLMyy(i,j3)
         End Do
        End Do

        Do i=1,Imax
         DLMxx(i,NSGQ) = DLMxx(i,NSGQ)-DLMxx(i,NSGT)
         DLMxy(i,NSGQ) = DLMxy(i,NSGQ)-DLMxy(i,NSGT)
         DLMyx(i,NSGQ) = DLMyx(i,NSGQ)-DLMyx(i,NSGT)
         DLMyy(i,NSGQ) = DLMyy(i,NSGQ)-DLMyy(i,NSGT)
        End Do

      end if

c--------------------
c Combine the incident flow
c with the slp
c--------------------

        Do i=1,Imax
         Uslp(i) = vsr1*(Uinf(i)+Uslp(i)/vs1)
         Vslp(i) = vsr1*(Vinf(i)+Vslp(i)/vs1)
c        write (6,100) i,Uslp(i),Uinf(i),Vslp(i),Vinf(i)
        End Do

c----------------------------
c testing integral identities
c----------------------------

      Ido = 1
      Ido = 0

      if(Ido.eq.1) then

      Do i=1,Imax

c        testxx = 0.0D0
c        testyx = 0.0D0
c        testxy = 0.0D0
c        testyy = 0.0D0

         testx  = 0.0D0    ! testing of eigenfunction
         testy  = 0.0D0    ! on a circular interface

         Do j=1,Imax

c         testxx = testxx+DLMxx(i,j)
c         testyx = testyx+DLMyx(i,j)
c         testxy = testxy+DLMxy(i,j)
c         testyy = testyy+DLMyy(i,j)

          uuuu = cos(TH2(j))   ! radial velocity
          vvvv = sin(TH2(j))   ! radial velocity
          testx = testx + DLMxx(i,j)*uuuu+ DLMyx(i,j)*vvvv
          testy = testy + DLMxy(i,j)*uuuu+ DLMyy(i,j)*vvvv

         End Do

c        write (6,100) i,x(i),y(i),testxx,testxy,testyx,testyy

         testx = cos(TH2(i)) - testx/pi2 + ttttx  ! should be zero
         testy = sin(TH2(i)) - testy/pi2 + tttty  ! should be zero

         if(Iflow.eq.6) then
           testx = testx + pi2*x(1) * x(i)
           testy = testy + pi2*x(1) * y(i)
         end if

         testm = sqrt(testx**2+testy**2)

         write (6,100) i,x(i),y(i),testx,testy,testm

        End Do

        end if

c-----------------------------
c prepare the iteration matrix
c-----------------------------

        Do i=1,Imax
         Do j=1,Imax

          DLMxx(i,j) = rkap2*DLMxx(i,j)
          DLMyx(i,j) = rkap2*DLMyx(i,j)
          DLMxy(i,j) = rkap2*DLMxy(i,j)
          DLMyy(i,j) = rkap2*DLMyy(i,j)

c         write (6,100) i,DLMxx(i,j),DLMyx(i,j),DLMxy(i,j),DLMyy(i,j)

         End Do
        End Do

c-------------------------------------
c THIS SECTION FOR AN EXPANDING BUBBLE
c-------------------------------------

      if(vs2.lt.0.000001) then

c-----
c generate the linear system
c for the node velocity
c-----


      Do i=1,Imax
       iImax = i+Imax
       Do j=1,Imax
         jImax = j+Imax
         Amat(i    ,j)      = -DLMxx(i,j)
         Amat(iImax,j)      = -DLMxy(i,j)
         Amat(i    ,j+Imax) = -DLMyx(i,j)
         Amat(iImax,j+Imax) = -DLMyy(i,j)
       End Do
       Amat(i,i)         = Amat(i,i)        +1.0D0
       Amat(iImax,iImax) = Amat(iImax,iImax)+1.0D0
       rhs(i)     = Uslp(i)
       rhs(iImax) = Vslp(i)
      End Do

c------------------
c Non-periodic flow
c------------------

      if(Iflow.ne.6) then

      call expand
     +
     +   (Imax
     +   ,Amat,rhs
     +   ,roexp
     +   ,mexp
     +   ,Iflow
     +   )

c--------------
c periodic flow
c--------------

      else 

       fc = roexp/tau

       Do i=1,Imax
         iImax = i+Imax
         rhs(i)     = rhs(i)     - fc*x(i)
         rhs(iImax) = rhs(iImax) - fc*y(i)
       End Do

c--- unecessary deflation

       if(mexp_foam.ne.0) then

       call expand
     +
     +   (Imax
     +   ,Amat,rhs
     +   ,roexp
     +   ,mexp
     +   ,Iflow
     +   )

       end if

c-----------
      end if
c-----------

c---
c solve the system
c---

      Mats = Imax2   ! matrix size

      Isym_g = 0     ! coefficient matrix is not symmetric
      Iwlpvt = 1     ! pivoting enabled

      call gel
     +
     +  (Mats
     +  ,Amat,rhs
     +  ,sln
     +  ,Isym_g
     +  ,Iwlpvt
c    +  ,l,u
     +  ,det
     +  ,Istop
     +  )

c-----------
c distribute the solution
c-----------

      Do i=1,Imax
       U(i) = sln(i)
       V(i) = sln(i+Imax)
c      write (6,100) i,U(i),V(i)
      End Do

      Go to 99

      End If

c-------------------------------------
c END OF SECTION FOR AN EXPANDING BUBBLE
c-------------------------------------

c----------------
c direct solution
c----------------

      if(Isolve.eq.1) then

      write (6,*) " caps_2d_vel: Solving the linear system "

      Do i=1,Imax

       IImax = i+Imax

       Do j=1,Imax
         Amat(i    ,j)      = -DLMxx(i,j)
         Amat(IImax,j)      = -DLMxy(i,j)
         Amat(i    ,j+Imax) = -DLMyx(i,j)
         Amat(IImax,j+Imax) = -DLMyy(i,j)
       End Do

       Amat(i,i)         = Amat(i,i)        +1.0D0
       Amat(IImax,IImax) = Amat(IImax,IImax)+1.0D0
       rhs(i)     = Uslp(i)
       rhs(IImax) = Vslp(i)

      End Do

      mats = 2*Imax  ! matrix size

      Isym_g = 0       ! coefficient matrix is not symmetric
      Iwlpvt = 1       ! pivoting enabled

      call gel
     +
     +  (mats
     +  ,Amat,rhs
     +  ,sln
     +  ,Isym_g
     +  ,Iwlpvt
c    +  ,l,u
     +  ,det
     +  ,Istop
     +  )

      Do i=1,Imax
       U(i) = sln(i)
       V(i) = sln(i+Imax)
      End Do

c---
c Done with the direct solution
c---

      write (6,*) " caps_2d_vel: linear system solved"

      Go to 99

      end if

c-------------------
c NEUMANN ITERATIONS
c-------------------

      write (6,*) "caps_2d_vel: Started Neumann Iterations"

      Iter = 0

  98  Continue

      Iter  = Iter + 1
      Iflag = 0

c     Do i=1,Imax
c       write (6,100) i,U(i),V(i)
c     End Do

c------------------------------------------------------
c
c  Deflation
c  ---------
c
c  Idfl = 0 no deflation
c  Idfl = 1 one eigenfunction will be removed 
c  Idfl = 2 all four eigenfunctions will be removed 
c
c-----------

                         ! four modes of deflation
      CF1 = 0.0          ! normal vector
      CF2 = 0.0          ! x translation
      CF3 = 0.0          ! y translation
      CF4 = 0.0          ! rotation about z

      if(Idfl.ne.0) then

c     write (6,*)
c     write (6,*) "Computing the deflation coefficients"
c     write (6,*)

        call deflation 
     +
     +   (Idfl
     +   ,Isym
     +   ,U,V
     +   ,CF1,CF2,CF3,CF4
     +   ,xcnt,ycnt
     +   )

       CF1 = rkap*CF1
       CF2 = rkap*CF2
       CF3 = rkap*CF3
       CF4 = rkap*CF4

      End If

c---------------------------------------------------------

c-----------
c  Updating
c-----------

      corr = 0.0D0

      Do 5 i=1,Imax

        Udl = 0.0D0
        Vdl = 0.0D0
        Do j=1,Imax
          Udl = Udl + DLMxx(i,j)*U(j) + DLMyx(i,j)*V(j)
          Vdl = Vdl + DLMxy(i,j)*U(j) + DLMyy(i,j)*V(j)
c         write (6,100) i,DLMxx(i,j),DLMyx(i,j)
c    +                   ,DLMxy(i,j),DLMyy(i,j)
        End Do

        Unew(i) = Uslp(i)+UDL 
        Vnew(i) = Vslp(i)+VDL 

c---
c deflation
c---

        Unew(i) = Unew(i)
     +          - CF1 * cos(th2(i))*ornt(i)
     +          + CF2
     +          + CF4*(-Y(i)+ycnt)

        Vnew(i) = Vnew(i)
     +          - CF1 * sin(th2(i))*ornt(i)
     +          + CF3
     +          + CF4*( X(i)-xcnt)

c---
c increments
c---

        DU = Unew(i)-U(i)
        DV = Vnew(i)-V(i)

        If(JGS.eq.2) then     ! Gauss-Siedel
          U(i) = Unew(i)
          V(i) = Vnew(i)
        End If

        If(abs(DU).ge.eps.or.abs(DV).ge.eps) Iflag = 1 
        If(abs(DU).ge.Corr) Corr = abs(DU)
        If(abs(DV).ge.Corr) Corr = abs(DV)

c       write (6,199) i,DU,DV,UDL,VDL

  5   Continue

c----------------
c End of updating
c----------------

c------
c diagnostics
c------

      if(Idfl.eq.0) write (6,195) Iter,Corr
      if(Idfl.eq.1) write (6,196) Iter,Corr,CF1
      if(Idfl.eq.2) write (6,197) Iter,Corr,CF1,CF2,CF3,CF4

c-------

      if(JGS.eq.1) then         ! Jacobi
        Do i=1,Imax
          U(i) = Unew(i)
          V(i) = Vnew(i)
        End Do
      end if

c-------------------------------------

      if(Iter.gt.Nter) then
        write (6,*)
        write (6,*) " Iterations failed to converge"
        Go to 91
      end if

      U(Nsg1) = U(1)     ! needed for deflation
      V(Nsg1) = V(1)     ! needed for deflation

      If(Iflag.eq.1) Go to 98 

c-------------------------
c Done with the iterations
c-------------------------

c----------------------------------------

  91  Continue

      If(Isolve.eq.1) Idfl = 0

c-------------------------------
c  if deflation of RBM was done,
c  rectify the solution
c-------------------------------

      if(Idfl.eq.2) then

       cf = 1.0/(1.0+rkap)
       Do i=1,Imax
         U(i) = U(i)-cf*(CF2
     +                  -CF4*(Y(i)-ycnt)
     +                  )
         V(i) = V(i)-cf*(CF3
     +                  +CF4*(X(i)-xcnt)
     +                  )
       End Do

      end if

c---------------------------------------

  99  Continue

c---
c reflect a symmetric interface
c---

      if(Isym.eq.2) then

        Do i=1,NSGQ
         U(NSGM-i+1) = -U(i)
         V(NSGM-i+1) =  V(i)
         U(NSGM+i-1) = -U(i)
         V(NSGM+i-1) = -V(i)
         U(NSG -i+2) =  U(i)
         V(NSG -i+2) = -V(i)
        End Do

      end if

c---
c wrap around
c---

      U(0) = U(NSG)
      V(0) = V(NSG)

      U(NSG1) = U(1)
      V(NSG1) = V(1)

      U(NSG2) = U(2)
      V(NSG2) = V(2)

c-----
c Done
c-----

 195  Format (' Iter:',I3,' Max Corr=',F9.6)
 196  Format (' Iter:',I3,' Max Corr=',F9.6,' CF1= ',F7.5)
 197  Format (' Iter:',I3,' Max Corr=',F9.6,' CF1= ',F7.5,
     +          ' CF2 = ',f7.5," CF3 =",F7.5," CF4=",F7.5)
 198  Format (' Evaluating the SLP at point:',I2)
 199  Format (' Increments at point ',I2,4(F10.5))
 100  Format (1X,I3,9(1X,F10.5))
 101  Format (1X,F15.10)

      return
      end
