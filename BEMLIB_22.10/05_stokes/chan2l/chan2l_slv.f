      subroutine chan2l_slv
     +
     +  (RL
     +  ,h
     +  ,NGL
     +  ,NSG
     +  ,rho1,vs1
     +  ,rho2,vs2
     +  ,V1,V2
     +  ,chi
     +  ,gx,gy
     +  ,vnx,vny
     +  ,crv
     +  ,srtn
     +  ,Ux,Uy
     +  ,Un,Ut
     +  )

c===========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

c-----------------------
c Two-layer channel flow
c
c Generate and solve
c a system of linear equations
c descending from the boundary-element
c collocation method
c
c NOTATION:
c -------
c
c RL:    period
c h:     channel semi-width
c NSG:   number of segments along the interface
c X, Y:  interfacial nodes
c-------------------------------------------------

      Implicit double precision (a-h,o-z)

      Dimension   X(0:513),  Y(0:513)
      Dimension vnx(0:513),vny(0:513)
      Dimension crv(0:513)

      Dimension DfxD(0:513),DfyD(0:513)

      Dimension   Ux(0:513),  Uy(0:513)
      Dimension   Un(0:513),  Ut(0:513)

      Dimension srtn(0:513)

      Dimension   X0(512),  Y0(512)
      Dimension DuR0(512)
      Dimension DfD0(512)

      Dimension  Ux0(512), Uy0(512)
      Dimension vnx0(512),vny0(512)

c---
c influence matrices
c---

      Dimension ASS(900,900),BSS(900,900),CM(900,900)
      Dimension rhs(900),sol(900)

c--------------
c common blocks
c--------------

      common/XXYY/X,Y

      common/ppii/pi,pih,pi2,pi4,pi6,pi8

c--------
c prepare
c--------

      NSG1 = NSG+1
      NSG2 = NSG+2

      hs = h**2
      xi = 0.50D0*(V2-V1)/h

      if(abs(V2-V1).gt.0.001) then
       yref = - h*(V2+V1)/(V2-V1)
      else
       yref = 0.0D0
      end if

      cf_slp = 1.0D0/(pi4*vs1)
      cf_dlp = 1.0D0/ pi4

      vsr = vs2/vs1
      vsr1 = vsr+1.0D0

      delta = rho2/rho1
      Drho  = rho1 - rho2

      cff = 2.0D0*(1.0D0-vsr)/vsr1

c--------------------------------------
c define the surface collocation points
c--------------------------------------

      Do i=1,NSG
        X0(i) = 0.50D0*( X(i)+X(i+1) )
        Y0(i) = 0.50D0*( Y(i)+Y(i+1) )
      End Do

c-----------------------------------------
c Generate the influence matrices ASS, BSS
c-----------------------------------------

      Do I=1,NSG    ! loop over surface collocation points

c     write (6,*)  " chan2l_slv: SS: collocating at wall point :",I

       NI = NSG+I

c---
c run over surface elements
c---

        Do J=1,NSG       ! loop over elements

         X1 = X(J)
         Y1 = Y(J)

         X2 = X(J+1)
         Y2 = Y(J+1)

         Ising = 0
         If(I.eq.J) Ising = 1   ! subtract off the singularity

         call chan2l_sdlp 
     +
     +      (RL
     +      ,h
     +      ,X0(I),Y0(I)
     +      ,X1,Y1
     +      ,X2,Y2
     +      ,NGL
     +      ,Ising
     +      ,Qxx,Qxy
     +      ,Qyx,Qyy
     +      ,Wxx,Wyx
     +      ,Wxy,Wyy
     +      )

         NJ = NSG +J

         ASS(I, J) = cf_slp*Qxx     ! x component of the BIE
         ASS(I,NJ) = cf_slp*Qyx     ! x component of the BIE

         ASS(NI, J) = cf_slp*Qxy     ! y component of the BIE
         ASS(NI,NJ) = cf_slp*Qyy     ! y component of the BIE

c-- minus sign because normal vector is downward

         BSS(I, J) = -cf_dlp*Wxx     ! x component of the BIE
         BSS(I,NJ) = -cf_dlp*Wyx     ! x component of the BIE

         BSS(NI, J) = -cf_dlp*Wxy     ! x component of the BIE
         BSS(NI,NJ) = -cf_dlp*Wyy     ! x component of the BIE

         End Do      ! run over elements

      End Do    ! run over collocation points

c--------------------------
c identities satisfied 
c by the influence matrices
c--------------------------

      Ido = 1
      Ido = 0

      if(Ido.eq.1) then

c---
c compute the normal vector
c---

      Do i=1,NSG
        i1 = i+1
        Dx = X(i1)-X(i)
        Dy = Y(i1)-Y(i)
        Dl = Dsqrt(Dx*Dx+Dy*Dy)
        vnx0(i) =  Dy/Dl
        vny0(i) = -Dx/Dl
      End Do

      Do i=1,2*NSG

       test1 = 0.0D0
       Do j=1,NSG
        test1 = test1 + ASS(i,j)    *vnx0(j)
     +                + ASS(i,j+NSG)*vny0(j)
       End Do

       test2 = 0.0D0
       Do j=1,NSG
        test2 = test2 + BSS(i,j)
c       write (6,*) i,j,BSS(i,j)
       End Do

       test3 = 0.0D0
       Do j=1,NSG
        test3 = test3 + BSS(i,NSG+j)
       End Do

       write (6,100) i,test1,test2,test3

      End Do

      pause

      end if

c------------------------------------
c Discontinuity of the reference flow
c at the element mid-points
c------------------------------------

      Do i=1,NSG
       DuR0(i) = 0.50D0/vs1 * ( chi * (1.0D0-1.0D0/vsr)
     +        + rho1*gx * (1.0D0-delta/vsr) ) *(hs-Y0(i)**2)
       DuR0(i+NSG) = 0.0D0
      End Do

c--------------------------------
c Compute Df at the surface nodes
c EXCLUDING the Marangoni stresses
c--------------------------------

c     write (6,*) "chan2l: vnx,vny,crv,srtn"

      Do i=1,NSG1

        vx = vnx(i)
        vy = vny(i)

        fc = Drho*Y(i)
        c1 = vs1*(1.0D0-vsr)*xi
        
        DfxR = - fc * (gy*vx + gx*vy) + c1 * vy
        DfyR = - fc * (gy*vy + gx*vx) + c1 * vx

        DfxD(i) = - srtn(i)*crv(i)*vx - DfxR
        DfyD(i) = - srtn(i)*crv(i)*vy - DfyR

c       write (6,100) i,vx,vy,crv(i),srtn(i)
c    +        ,DfxD(i),DfyD(i),DfxR,DfyR
      End Do

c     pause

c---------------------------------
c compute Df and normal vector
c at the surface collocation points
c
c Include the Marangoni tractions
c---------------------------------

      Do i=1,NSG

        i1 = i+1

        Dx = X(i1)-X(i)
        Dy = Y(i1)-Y(i)
        Dl = Dsqrt(Dx*Dx+Dy*Dy)
        Dsrtn = srtn(i1)-srtn(i)
        tx = Dx/Dl
        ty = Dy/Dl

        DsrtnDl = Dsrtn/Dl

        DfD0(i)     = 0.50D0*( DfxD(i)+DfxD(i1) ) - DsrtnDl * tx
        DfD0(i+NSG) = 0.50D0*( DfyD(i)+DfyD(i1) ) - DsrtnDl * ty

c       write (6,*) i,X0(i),DsrtnDl
c       write (6,100) i,DfD0(i),DfD0(i+NSG)

      End Do

c     pause

c----------------------------------
c Generate the master linear system
c
c   AMLS . q = rhs = BMLS . f
c
c the unknown vector consists of the:
c
c  wall traction
c surface velocity
c----------------------------------

      Nsys = 2*NSG

c     write (6,*) " chan2l_slv: generating the MLS"

c--------------------------------
c generate the coefficient matrix
c--------------------------------

      Do i=1,Nsys
       Do j=1,Nsys
        CM(i,j) = - cff*BSS(i,j)
       End Do
        CM(i,i) = CM(i,i) + 1.0D0
      End Do

c-----------------
c generate the RHS
c-----------------

      Do i=1,Nsys

        sum = 0.0D0
        Do j=1,Nsys
         sum = sum + BSS(i,j)*DuR0(j)
        End Do

        rhs(i) = - vsr* (DuR0(i)+ 2.0D0 * sum)/vsr1

        sum = 0.0D0
        Do j=1,Nsys
         sum = sum + ASS(i,j)*DfD0(j)
        End Do

        rhs(i) = rhs(i) - 2.0D0*sum/vsr1

      End Do 

c     write (6,*) " chan2l_slv: MLS assembled"

c--------------
c solve the MLS
c--------------

c     write (6,*) " chan2l_slv: Solving the MLS"

      Isym_gl = 0
      Iwlpvt_gl = 1

      call gel
     +
     +   (Nsys,CM
     +   ,RHS
     +   ,sol
     +   ,Isym_gl
     +   ,Iwlpvt_gl
c    +   ,l,u
     +   ,det
     +   ,Istop
     +   )

c---------
c printing
c---------

c     Do i=1,Nsys
c       write (6,100) i,rhs(i)
c       write (6,100) i,sol(i)
c       write (6,102) (CM(i,j),j=1,Nsys),RHS(i),sol(i)
c     End Do

c------------------------
c distribute the solution
c------------------------

      Ic = 0

      Do i=1,NSG
       Ic = Ic+1
       Ux0(i) = sol(Ic)
c      Unp =    xi *(Y0(i)-yref) 
c    +      + 0.5D0*(chi+rho1*gx)/vs1 * (hs-Y0(i)**2)
c      Ux0(i) = Ux0(i)+Unp
c      write (6,100) Ic,sol(Ic),Unp
      End Do

      Do i=1,NSG 
       Ic = Ic+1
       Uy0(i) = sol(Ic)
c      write (6,*) Ic,sol(Ic)
      End Do

c--------------------------
c compute the node velocity
c--------------------------

      Ux(1) = 0.50D0*(Ux0(NSG)+Ux0(1))
      Uy(1) = 0.50D0*(Uy0(NSG)+Uy0(1))

      Do i=2,NSG
       Ux(i) = 0.50D0*( Ux0(i-1)+Ux0(i) )
       Uy(i) = 0.50D0*( Uy0(i-1)+Uy0(i) )
      End Do

      Do i=1,NSG
        Unp =  xi *(Y(i)-yref) 
     +      + 0.5D0*(chi+rho1*gx)/vs1 * (hs-Y(i)*Y(i))
       Ux(i) = Ux(i)+Unp
c      write (6,100) i,Ux(i),Uy(i),Unp
      End Do

c     stop

c----------------------------------
c compute the
c normal and tangential components
c of the interfacial node velocity
c----------------------------------

      Do i=1,NSG
       tx = -vny(i)
       ty =  vnx(i)
       Un(i) = Ux(i)*vnx(i) + Uy(i)*vny(i)
       Ut(i) = Ux(i)*tx     + Uy(i)*ty
c      write (6,103) i,Ux(i),Uy(i),Un(i),Ut(i)
      End Do

      Ux(0) = Ux(NSG)
      Uy(0) = Uy(NSG)
      Un(0) = Un(NSG)
      Ut(0) = Ut(NSG)

      Ux(NSG1) = Ux(1)
      Uy(NSG1) = Uy(1)
      Un(NSG1) = Un(1)
      Ut(NSG1) = Ut(1)

      Ux(NSG2) = Ux(2)
      Uy(NSG2) = Uy(2)
      Un(NSG2) = Un(2)
      Ut(NSG2) = Ut(2)

c-----
c done
c-----

c     write (6,*) "chan2l_slv: exiting"

  100 Format (1X,I4,10(1x,F10.7))
  101 Format (10(1x,F15.10))
  102 Format (200(1x,F6.4))
  103 Format (1X,I4,10(1x,F10.5))
  109 Format (200(1x,F10.3))

      return
      end
