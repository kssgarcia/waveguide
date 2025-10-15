      subroutine surfact_fem
     +
     +   (Nelm
     +   ,Npts
     +   ,mint
     +   ,Dt
     +   ,Ds
     +   ,Move
     +   )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c---------------------------------
c compute the global:
c
c  conductivity matrix
c  advection matrix
c  mass matrix
c
c
c LEGEND:
c ------

c edf (6,6): element diffusion matrix
c ead (6,6): element advection matrix
c ema (6,6): element mass matrix
c emd (6,6): element mass matrix weighted by theta
c emk (6,6): element mass matrix weighted by 2 kappa u_n
c
c ph(i): element interpolation function psi, i = 1, ..., 6
c
c gph(i,3): surface gradient of psi(i),
c---------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension      p(1026,3)
      Dimension     ne(1026,7)
      Dimension    vna(1026,3)
      Dimension      u(1026,3)
      Dimension      c(1026)

      Dimension     n(512,6),nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension crvm(1026)
      Dimension dilt(1026)

      Dimension xiq(20),etq(20),wq(20)

      Dimension  ph(6),  gph(6,3)
      Dimension edf(6,6),ead(6,6),ema(6,6)
      Dimension emd(6,6),emk(6,6)

      Dimension Gdf(1026,1026),Gad(1026,1026),Gma(1026,1026)
      Dimension Gmd(1026,1026),Gmk(1026,1026)

      Dimension  B(1026,1026)
      Dimension AM(1026,1026),rhs(1026)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/geo2/vna
      common/geo3/crvm

      common/surfa/c

      common/veloc0/u

      common/trq/xiq,etq,wq

c--------
c prepare
c--------

      Dth = 0.50D0*Dt

c------------------------------------------------
c When Move = 0 (points move with fluid velocity)
c
c compute the rate of dilatation
c due to in-plane stretching :
c
c grad_s . u_s
c
c-----------

      if(Move.eq.0) then

c       write (6,*) " surfact_fvm: dilatation due to in-plane"
c       write (6,*) "              stretching"

        call dilt_3d
     +
     +  (Nelm
     +  ,Npts
     +  ,Dilt
     +  )

c       Do i=1,Npts
c        write (6,*) i,Dilt(i)
c       End Do
c       pause

      end if

c-------------------------------
c initialize the global matrices
c------------------------------ 

      Do k=1,Npts
        Do l=1,Npts
          Gdf(k,l) = 0.0D0     ! diffusion
          Gad(k,l) = 0.0D0
          Gma(k,l) = 0.0D0
          Gmd(k,l) = 0.0D0
          Gmk(k,l) = 0.0D0
        End Do
      End Do

c-----------------------
c loop over the elements
c-----------------------

      Do k=1,Nelm

      i1 = n(k,1)   ! global element node labels
      i2 = n(k,2)
      i3 = n(k,3)
      i4 = n(k,4)
      i5 = n(k,5)
      i6 = n(k,6)

      al = alpha(k)
      be = beta (k)
      ga = gamma(k)

c----------------------------------
c launch the integration quadrature
c----------------------------------

       Do j1=1,6
        Do j2=1,6
          edf(j1,j2) = 0.0D0
          ead(j1,j2) = 0.0D0
          ema(j1,j2) = 0.0D0
          emd(j1,j2) = 0.0D0
          emk(j1,j2) = 0.0D0
        End Do
       End Do

      Do i=1,mint

        xi  = xiq(i)
        eta = etq(i)

        call surfact_fem_interp
     +
     +  (p(i1,1),p(i1,2),p(i1,3)
     +  ,p(i2,1),p(i2,2),p(i2,3)
     +  ,p(i3,1),p(i3,2),p(i3,3)
     +  ,p(i4,1),p(i4,2),p(i4,3)
     +  ,p(i5,1),p(i5,2),p(i5,3)
     +  ,p(i6,1),p(i6,2),p(i6,3)
     +
     +  ,dilt(i1),dilt(i2),dilt(i3)
     +  ,dilt(i4),dilt(i5),dilt(i6)
     +
     +  ,crvm(i1),crvm(i2),crvm(i3)
     +  ,crvm(i4),crvm(i5),crvm(i6)
     +
     +  ,vna(i1,1),vna(i1,2),vna(i1,3)
     +  ,vna(i2,1),vna(i2,2),vna(i2,3)
     +  ,vna(i3,1),vna(i3,2),vna(i3,3)
     +  ,vna(i4,1),vna(i4,2),vna(i4,3)
     +  ,vna(i5,1),vna(i5,2),vna(i5,3)
     +  ,vna(i6,1),vna(i6,2),vna(i6,3)
     +
     +  ,u(i1,1),u(i1,2),u(i1,3)
     +  ,u(i2,1),u(i2,2),u(i2,3)
     +  ,u(i3,1),u(i3,2),u(i3,3)
     +  ,u(i4,1),u(i4,2),u(i4,3)
     +  ,u(i5,1),u(i5,2),u(i5,3)
     +  ,u(i6,1),u(i6,2),u(i6,3)
     +
     +  ,al,be,ga
     +  ,xi,eta
     +
     +  ,diltp
     +  ,crvmp
     +  ,vnx,vny,vnz
     +  ,ux,uy,uz
     +  ,ph,gph
     +  ,hs
     +  )

       cf = 0.5D0*hs*wq(i)

       Do j1=1,6
        Do j2=1,6
c       Do j2=j1,6

        prj =  gph(j1,1)*gph(j2,1)   
     +       + gph(j1,2)*gph(j2,2)
     +       + gph(j1,3)*gph(j2,3)

        edf(j1,j2) = edf(j1,j2) + prj * cf          ! diffusion

        ema(j1,j2) = ema(j1,j2) + ph(j1)*ph(j2)*cf  ! mass 

c-------------------------
        if(Move.eq.0) then
c-------------------------

         emd(j1,j2) = emd(j1,j2) + ph(j1)*ph(j2)*diltp*cf

c------------------------------
        else if(Move.eq.1) then
c------------------------------

         ead(j1,j2) = ead(j1,j2)
     +       + ( ux * gph(j1,1) + uy * gph(j1,2)
     +          +uz * gph(j1,3) )*ph(j2)*cf      ! advection

         emk(j1,j2) = emk(j1,j2)
     +        + ph(j1)*ph(j2)*2.0D0*crvmp
     +          *(ux*vnx + uy*vny + uz*vnz)*cf
c-------------
        end if
c-------------

        End Do
       End Do

      End Do

c-----
c symmetric components
c------
c
c     Do j1=2,6
c      Do j2=1,j1-1
c        edf(j1,j2) = edf(j2,j1)
c        ema(j1,j2) = ema(j2,j1)
c      End Do
c     End Do

c------------------------------------------
c quadrature finished
c
c make contributions to the global matrices
c------------------------------------------

      Do k1=1,6

       if(k1.eq.1) l1 = i1
       if(k1.eq.2) l1 = i2
       if(k1.eq.3) l1 = i3
       if(k1.eq.4) l1 = i4
       if(k1.eq.5) l1 = i5
       if(k1.eq.6) l1 = i6

       Do k2=1,6

        if(k2.eq.1) l2 = i1
        if(k2.eq.2) l2 = i2
        if(k2.eq.3) l2 = i3
        if(k2.eq.4) l2 = i4
        if(k2.eq.5) l2 = i5
        if(k2.eq.6) l2 = i6
 
        Gdf(l1,l2) = Gdf(l1,l2) + edf(k1,k2)
        Gma(l1,l2) = Gma(l1,l2) + ema(k1,k2)

        if(Move.eq.0) then
          Gmd(l1,l2) = Gmd(l1,l2) + emd(k1,k2)
        else if(Move.eq.1) then
          Gad(l1,l2) = Gad(l1,l2) + ead(k1,k2)
          Gmk(l1,l2) = Gmk(l1,l2) + emk(k1,k2)
        end if

       End Do
      End Do

c------------------------------
c End of loop over the elements
c------------------------------

      End Do

c--------------------
c Define the B matrix
c--------------------

c----------------------
      if(Move.eq.0) then   ! total velocity
c----------------------

      Do i=1,Npts
       Do j=1,Npts
        B(i,j) = -Gmd(i,j) - Ds * Gdf(i,j)
c       write (6,102) i,j,B(i,j)
       End Do
      End Do

c----------------------------
      else if(Move.eq.1) then  ! normal velocity
c----------------------------

      Do i=1,Npts
       Do j=1,Npts
        B(i,j) = Gad(i,j) - Gmk(i,j) - Ds * Gdf(i,j)
c       write (6,102) i,j,B(i,j)
       End Do
      End Do

c-----------
      end if
c-----------

c---------------------------
c complete the linear system
c---------------------------

      Do i=1,Npts
       Do j=1,Npts
        B(i,j)=Dth*B(i,j)
       End Do
      End Do

c------------------------------------------
c compute the RHS of the transport equation
c------------------------------------------

      Do i=1,Npts
       RHS(i) = 0.0D0
       Do j=1,Npts
        RHS(i) = RHS(i) + (Gma(i,j)+B(i,j))*c(j)
       End Do
      End Do

c---------------
c compute the coefficient matrix
c of the transport equation
c---------------

      Do i=1,Npts
       Do j=1,Npts
        AM(i,j) = Gma(i,j) - B(i,j)
       End Do
      End Do

c------------------------
c Solve the linear system
c------------------------

      write (6,*) "surfact_fem: Solving the linear system ", Npts

      Isym_gel   = 0
      Iwlpvt_gel = 1

      call gel
     +
     + (Npts
     + ,AM
     + ,RHS
     + ,c
     + ,Isym_gel
     + ,Iwlpvt_gel
c    + ,l,u
c    + ,det
     + ,Istop
     + )

c-----
c done
c-----

 102  Format (1x,2(1x,i3),2(1x,f15.10))

      return
      end
