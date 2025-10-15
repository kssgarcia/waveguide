      subroutine ldlp_3d_se_integral
     +
     +     (x0,y0,z0
     +     ,nsp
     +     ,psp
     +     ,k
     +     ,intm
     +     ,mint
     +     ,q
     +     ,q0
     +     ,ptl
     +     ,area
     +     ,qsint
     +     )

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c--------------------------------------
c Compute the double-layer laplace potential
c over a non-singular triangle numbered k
c--------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension   p(1026,3)
      Dimension  ne(1026,7)
      Dimension vna(1026,3)

      Dimension     n(512,6),nbe(512,3)
      Dimension alpha(512),beta(512),gamma(512)

      Dimension xiq(20),etq(20),wq(20)

c spectral:

      Dimension vmaster(8),van(500,500),vaninv(500,500)
      Dimension   xsp(512,100),ysp(512,100),zsp(512,100)

      Dimension q(2000)
      Dimension   nsp(512,100)   ! connectivity
      Dimension   psp(20000,3)

      Dimension aux(100),psisp(100)

      Dimension xir(250000),etr(250000),wr(250000)    ! integration weights

      common/spectral1/mpoly,Npoly,vmaster,van,vaninv
      common/spectral2/nesp,ngsp

c---
c common blocks
c---

      common/points/p,ne
      common/elmnts/n,nbe

      common/albega/alpha,beta,gamma

      common/geo2/vna

      common/var1/Iflow,Ign,wall

      common/trq/xiq,etq,wq

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c----------
c set flags
c----------

      Iopt_int = 2   ! need the surface metric
      Iopt_lgf = 2   ! need grad(G)

c-----------
c initialize
c-----------

      area  = 0.0D0
      qsint = 0.0D0
      ptl   = 0.0D0

c---
c define triangle nodes
c and mapping parameters
c---

      i1 = n(k,1)
      i2 = n(k,2)
      i3 = n(k,3)
      i4 = n(k,4)
      i5 = n(k,5)
      i6 = n(k,6)

      al = alpha(k)
      be = beta (k)
      ga = gamma(k)

c-----------------------------
c loop over integration points
c-----------------------------

      if(intm.eq.1) then    ! first integration method: quadrature

        mintt = mint

      else if(intm.eq.2) then   ! second integration method: Newton-Cotes

       call trgl_rule (mint,Nbase,xir,etr,wr)

c      write (6,*) Nbase
c      write (6,*)
c      write (6,*) (xir(i),i=1,Nbase)
c      write (6,*)
c      write (6,*) (etr(i),i=1,Nbase)
c      write (6,*)
c      write (6,*) ( wr(i),i=1,Nbase)
c      stop

       mintt = Nbase

      end if

c---
c base points and weights
c---

      Do 1 i=1,mintt    ! quadrature

        if(intm.eq.1) then
          xi  = xiq(i)    ! quadrature
          eta = etq(i)     ! quadrature
        else if(intm.eq.2) then
          xi  = xir(i)    ! quadrature
          eta = etr(i)     ! quadrature
        end if

        if(xi.lt.0.000001)  xi=0.001D0
        if(xi.gt.0.999999)  xi=0.999D0
        if(eta.lt.0.000001) eta=0.001D0
        if(eta.gt.0.999999) eta=0.999D0

        call ldlp_3d_se_interp
     +
     +    (Iopt_int
     +    ,p(i1,1),p(i1,2),p(i1,3)
     +    ,p(i2,1),p(i2,2),p(i2,3)
     +    ,p(i3,1),p(i3,2),p(i3,3)
     +    ,p(i4,1),p(i4,2),p(i4,3)
     +    ,p(i5,1),p(i5,2),p(i5,3)
     +    ,p(i6,1),p(i6,2),p(i6,3)
     +
     +    ,vna(i1,1),vna(i1,2),vna(i1,3)
     +    ,vna(i2,1),vna(i2,2),vna(i2,3)
     +    ,vna(i3,1),vna(i3,2),vna(i3,3)
     +    ,vna(i4,1),vna(i4,2),vna(i4,3)
     +    ,vna(i5,1),vna(i5,2),vna(i5,3)
     +    ,vna(i6,1),vna(i6,2),vna(i6,3)
     +
     +     ,al,be,ga
     +     ,xi,eta
     +
     +     ,x,y,z
     +     ,vnx,vny,vnz
     +     ,hs
     +     )

c---
c compute the Green's function
c---

c------
      if(Iflow.eq.1) then
c------

       call lgf_3d_fs 
     +
     +    (Iopt_lgf
     +    ,x,y,z
     +    ,x0,y0,z0
     +    ,G
     +    ,Gx,Gy,Gz
     +    )

c------
      else if(Iflow.eq.2) then
c------

       call lgf_3d_w
     +
     +   (Iopt_lgf
     +   ,Ign
     +   ,x,y,z
     +   ,x0,y0,z0
     +   ,wall
     +   ,G
     +   ,Gx,Gy,Gz
     +   )

c------
       end if
c------

c        If(abs(G).gt.100000.0) then
c        write (6,100) x,y,z,x0,y0,z0,G,xi,eta
c        End If

c-----
c interpolate for the dipole
c using the Vandermonde matrix:
c-----

         xip = 2.0D0*xi/(1.D0-eta)-1.0D0
         etap = 2.0D0*eta-1.0D0

         Jc=0
         Do ksp=0,mpoly
          Do lsp=0,mpoly-ksp
            Jc=Jc+1
            call jacobi (0D0,0D0,ksp,xip,rjac1)
            call jacobi (2.0D0*ksp+1.0D0,0.0D0,lsp,etap,rjac2)
            aux(Jc) = rjac1 *(1.0D0-eta)**ksp *rjac2
c           write (6,*) aux(Jc)
           End Do
         End Do

        Do isp=1,nesp
         psisp(isp) = 0.0D0
          Do jsp=1,nesp
            psisp(isp) = psisp(isp) + vaninv(isp,jsp)*aux(jsp)
c           write (6,*) vaninv(jsp,isp)
          End Do
        End Do

c   construct q by Lagrange intepolation:

        qint = 0.0D0

        Do isp=1,nesp
         qint = qint + psisp(isp)*q(nsp(k,isp))
c        write (6,*) psisp(isp),q(nsp(k,isp))
        End Do

c---
c apply the triangle quadrature
c---

      if(intm.eq.1) then        ! quadrature
       cf = 0.5D0*hs*wq(i)   
      else if(intm.eq.2) then
       cf = 0.5D0*hs*wr(i)     ! integration rule
      end if

      area = area + cf
      qsint = qsint + qint*cf

      ptl = ptl + (qint-q0)*(vnx*Gx+vny*Gy+vnz*Gz)*cf

c     write (6,100) G,xi,eta,x,y,z,x0,y0,z0


  1   Continue    ! end of loop over integration points

c-----
c done
c-----

 100  Format (100(1x,f10.5))

      return
      end
