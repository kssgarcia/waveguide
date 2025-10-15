      subroutine surfact_fvm
     +
     +  (Nelm   ! number of elements
     +  ,Npts   ! number of points
     +  ,mint
     +  ,Dt     ! time step
     +  ,Ds
     +  ,cel    ! element concentration
     +  ,Move
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

c------------------------------------------
c Advance the surfactant concentration using
c an implicit method
c
c Move = 0:  points move with total velocity
c        1:  points move with normal velocity
c------------------------------------------

      Implicit double precision (a-h,o-z)

      Dimension     p(1026,3)
      Dimension    ne(1026,7)
      Dimension   vna(1026,3)
      Dimension     u(1026,3)
      Dimension     c(1026)
      Dimension  dilt(1026)

      Dimension crvm(1026)

      Dimension      n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)
      Dimension arel(512)
      Dimension xmom(512),ymom(512),zmom(512)
      Dimension  cel(512)

      Dimension Dexpel(512)
      Dimension Diltel(512)

      Dimension  c_imp(1026),cel_imp(512)
      Dimension  gradc_imp(1026,3)

      Dimension  xxi(6), eet(6)
      Dimension DxDx(6),DyDx(6),DzDx(6)
      Dimension DxDe(6),DyDe(6),DzDe(6)
      Dimension   vx(6),  vy(6),  vz(6)
      Dimension rntx(6),rnty(6),rntz(6)

      Dimension xiq(20),etq(20),wq(20)

      Dimension  A(512,512)

      Dimension AM(1026,1026),rhs(1026),sol(1026)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/geo1/arel
      common/geo2/vna
      common/geo3/crvm
      common/geo9/xmom,ymom,zmom

      common/surfa/c

      common/veloc0/u
      common/trq/xiq,etq,wq

c---------
c prepare
c--------

      write (6,*) " surfact_fvm: Move =",Move
      write (6,*) " surfact_fvm: Dt   =",Dt

c----
      if(Move.eq.1) then  ! normal motion
c----

c---------------------------------------
c Nodes move normal to the interface
c compute the element rate of dilatation
c due to expansion: Int[ crvm u.n dS]
c---------------------------------------

c     write (6,*) " surfact_fvm: surface dilatation"
c     write (6,*) "              due to expansion"

      Do k=1,Nelm    ! loop over elements

      i1 = n(k,1)
      i2 = n(k,2)
      i3 = n(k,3)
      i4 = n(k,4)
      i5 = n(k,5)
      i6 = n(k,6)

      al = alpha(k)
      be =  beta(k)
      ga = gamma(k)

      Dexpel(k) = 0.0D0

      Do i=1,mint

        xi  = xiq(i)
        eta = etq(i)

c------------------------------
c Interpolate for the velocity, normal vector,
c and curvature
c at triangle quadrature points
c------------------------------

        call surfact_fvm_interp1
     +
     +  (p(i1,1),p(i1,2),p(i1,3)
     +  ,p(i2,1),p(i2,2),p(i2,3)
     +  ,p(i3,1),p(i3,2),p(i3,3)
     +  ,p(i4,1),p(i4,2),p(i4,3)
     +  ,p(i5,1),p(i5,2),p(i5,3)
     +  ,p(i6,1),p(i6,2),p(i6,3)
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
     +  ,x,y,z
     +  ,crvmp
     +  ,vnx,vny,vnz
     +  ,ux,uy,uz
     +  ,hs
     +  )

       cf = 0.5D0*hs*wq(i)

       udn = ux*vnx+uy*vny+uz*vnz

       Dexpel(k) = Dexpel(k)+ crvmp * udn * cf

      End Do

      Dexpel(k) = 2.0D0*Dexpel(k)/arel(k)

c     write (6,100) k,Dexpel(k)

      End Do

c----
      elseif(Move.eq.0) then
c----

c---------------------------------------
c Nodes move with the fluid velocity
c compute the element rate of dilatation
c due to in-plane stretching 
c
c grad_s . u_s
c
c Only when Move = 0
c---------------------------------------

c     write (6,*) " surfact_fvm: dilatation due to in-plane"
c     write (6,*) "              stretching"

      call dilt_3d
     +
     +  (Nelm
     +  ,Npts
     +  ,Dilt
     +  )

      Do k=1,Nelm

       Diltel(k) = 0.0D0

       Do j=1,6
        i = n(k,j)
        Diltel(k) = Diltel(k) + Dilt(i)
       End Do

       Diltel(k) = Diltel(k)/6.0D0

      End Do

c----
      end if
c----

c-------------------------------------
c  initialize temporary concentrations
c-------------------------------------

      Do i=1,Npts        ! nodes
       c_imp(i) = 0.0D0
      End Do

      Do i=1,Nelm        ! elements
       cel_imp(i) = 0.0D0
      End Do

c--------------------
c--------------------
c  begin the impulses
c--------------------
c--------------------

      Do 10 l=1,Nelm      ! IMPULSES

       cel_imp(l) = 1.0D0

c-------------------------
c Distribute concentration
c from elements to nodes
c-------------------------

      call surfact_fvm_etn
     +
     +  (Npts
     +  ,cel_imp
     +  ,c_imp
     +  )

c----------------------------------
c compute the surface grad of c_imp
c----------------------------------

      call sgrad_3d
     +
     + (Npts,Nelm
     + ,c_imp
     + ,gradc_imp
     + )

c-------------------------------------
c run over elements
c compute the element contour integral
c-------------------------------------

      Do 20 k=1,Nelm          ! RUN OVER THE ELEMENTS

       i1 = n(k,1)
       i2 = n(k,2)
       i3 = n(k,3)
       i4 = n(k,4)
       i5 = n(k,5)
       i6 = n(k,6)

       al = alpha(k)
       be =  beta(k)
       ga = gamma(k)

       alc = 1.0D0-al
       bec = 1.0D0-be
       gac = 1.0D0-ga

       xxi(1) = 0.0D0
       eet(1) = 0.0D0
       xxi(2) = 1.0D0
       eet(2) = 0.0D0
       xxi(3) = 0.0D0
       eet(3) = 1.0D0
       xxi(4) = al
       eet(4) = 0.0D0
       xxi(5) = ga
       eet(5) = gac
       xxi(6) = 0.0D0
       eet(6) = be

c----
c interpolate for 6 node values
c----

       Do i=1,6

        xi  = xxi(i)
        eta = eet(i)

        call surfact_fvm_interp
     +
     +  (p(i1,1),p(i1,2),p(i1,3)
     +  ,p(i2,1),p(i2,2),p(i2,3)
     +  ,p(i3,1),p(i3,2),p(i3,3)
     +  ,p(i4,1),p(i4,2),p(i4,3)
     +  ,p(i5,1),p(i5,2),p(i5,3)
     +  ,p(i6,1),p(i6,2),p(i6,3)
     +
     +  ,u(i1,1),u(i1,2),u(i1,3)
     +  ,u(i2,1),u(i2,2),u(i2,3)
     +  ,u(i3,1),u(i3,2),u(i3,3)
     +  ,u(i4,1),u(i4,2),u(i4,3)
     +  ,u(i5,1),u(i5,2),u(i5,3)
     +  ,u(i6,1),u(i6,2),u(i6,3)
     +
     +  ,c_imp(i1),c_imp(i2),c_imp(i3)
     +  ,c_imp(i4),c_imp(i5),c_imp(i6)
     +
     +  ,gradc_imp(i1,1),gradc_imp(i1,2),gradc_imp(i1,3)
     +  ,gradc_imp(i2,1),gradc_imp(i2,2),gradc_imp(i2,3)
     +  ,gradc_imp(i3,1),gradc_imp(i3,2),gradc_imp(i3,3)
     +  ,gradc_imp(i4,1),gradc_imp(i4,2),gradc_imp(i4,3)
     +  ,gradc_imp(i5,1),gradc_imp(i5,2),gradc_imp(i5,3)
     +  ,gradc_imp(i6,1),gradc_imp(i6,2),gradc_imp(i6,3)
     +
     +  ,al,be,ga
     +  ,xi,eta
     +  
     +  ,Move  
     +
     +  ,x,y,z
     +  ,ux,uy,uz
     +  ,cint
     +  ,gcx,gcy,gcz
     +  ,DxDx(i),DyDx(i),DzDx(i)
     +  ,DxDe(i),DyDe(i),DzDe(i)
     +  ,vx(i),vy(i),vz(i)
     +  ,hxi,het,hs
     +  )

       rntx(i) = Ds*gcx
       rnty(i) = Ds*gcy
       rntz(i) = Ds*gcz

       if(Move.eq.1) then
        rntx(i) = rntx(i) - cint * ux
        rnty(i) = rnty(i) - cint * uy
        rntz(i) = rntz(i) - cint * uz
       end if

       End Do

c-----------------------------------
c line integral along segment 1-4-2 (eta = 0)
c-----------------------------------

      bvx1 = vy(1)*DzDx(1)-vz(1)*DyDx(1)
      bvy1 = vz(1)*DxDx(1)-vx(1)*DzDx(1)
      bvz1 = vx(1)*DyDx(1)-vy(1)*DxDx(1)

      bvx2 = vy(4)*DzDx(4)-vz(4)*DyDx(4)
      bvy2 = vz(4)*DxDx(4)-vx(4)*DzDx(4)
      bvz2 = vx(4)*DyDx(4)-vy(4)*DxDx(4)

      bvx3 = vy(2)*DzDx(2)-vz(2)*DyDx(2)
      bvy3 = vz(2)*DxDx(2)-vx(2)*DzDx(2)
      bvz3 = vx(2)*DyDx(2)-vy(2)*DxDx(2)

      prj1 = bvx1*rntx(1) + bvy1*rnty(1) + bvz1*rntz(1)
      prj2 = bvx2*rntx(4) + bvy2*rnty(4) + bvz2*rntz(4)
      prj3 = bvx3*rntx(2) + bvy3*rnty(2) + bvz3*rntz(2)

      collect = al*prj1 + prj2 + alc*prj3

c-----------------------------------
c line integral along segment 2-5-3 (mixed)
c-----------------------------------

      bvx1 = vy(2)*DzDx(2)-vz(2)*DyDx(2)
      bvy1 = vz(2)*DxDx(2)-vx(2)*DzDx(2)
      bvz1 = vx(2)*DyDx(2)-vy(2)*DxDx(2)

      bvx2 = vy(5)*DzDx(5)-vz(5)*DyDx(5)
      bvy2 = vz(5)*DxDx(5)-vx(5)*DzDx(5)
      bvz2 = vx(5)*DyDx(5)-vy(5)*DxDx(5)

      bvx3 = vy(3)*DzDx(3)-vz(3)*DyDx(3)
      bvy3 = vz(3)*DxDx(3)-vx(3)*DzDx(3)
      bvz3 = vx(3)*DyDx(3)-vy(3)*DxDx(3)

      prj1 = bvx1*rntx(2) + bvy1*rnty(2) + bvz1*rntz(2)
      prj2 = bvx2*rntx(5) + bvy2*rnty(5) + bvz2*rntz(5)
      prj3 = bvx3*rntx(3) + bvy3*rnty(3) + bvz3*rntz(3)

      collect = collect - gac*prj1 - prj2 - ga*prj3

      bvx1 = vy(2)*DzDe(2)-vz(2)*DyDe(2)
      bvy1 = vz(2)*DxDe(2)-vx(2)*DzDe(2)
      bvz1 = vx(2)*DyDe(2)-vy(2)*DxDe(2)

      bvx2 = vy(5)*DzDe(5)-vz(5)*DyDe(5)
      bvy2 = vz(5)*DxDe(5)-vx(5)*DzDe(5)
      bvz2 = vx(5)*DyDe(5)-vy(5)*DxDe(5)

      bvx3 = vy(3)*DzDe(3)-vz(3)*DyDe(3)
      bvy3 = vz(3)*DxDe(3)-vx(3)*DzDe(3)
      bvz3 = vx(3)*DyDe(3)-vy(3)*DxDe(3)

      prj1 = bvx1*rntx(2) + bvy1*rnty(2) + bvz1*rntz(2)
      prj2 = bvx2*rntx(5) + bvy2*rnty(5) + bvz2*rntz(5)
      prj3 = bvx3*rntx(3) + bvy3*rnty(3) + bvz3*rntz(3)

      collect = collect + gac*prj1 + prj2 + ga*prj3

c------------------------------------------
c line integral along segment 3-6-1 (xi = 0)
c------------------------------------------

      bvx1 = vy(1)*DzDe(1)-vz(1)*DyDe(1)
      bvy1 = vz(1)*DxDe(1)-vx(1)*DzDe(1)
      bvz1 = vx(1)*DyDe(1)-vy(1)*DxDe(1)

      bvx2 = vy(6)*DzDe(6)-vz(6)*DyDe(6)
      bvy2 = vz(6)*DxDe(6)-vx(6)*DzDe(6)
      bvz2 = vx(6)*DyDe(6)-vy(6)*DxDe(6)

      bvx3 = vy(3)*DzDe(3)-vz(3)*DyDe(3)
      bvy3 = vz(3)*DxDe(3)-vx(3)*DzDe(3)
      bvz3 = vx(3)*DyDe(3)-vy(3)*DxDe(3)

      prj1 = bvx1*rntx(1) + bvy1*rnty(1) + bvz1*rntz(1)
      prj2 = bvx2*rntx(6) + bvy2*rnty(6) + bvz2*rntz(6)
      prj3 = bvx3*rntx(3) + bvy3*rnty(3) + bvz3*rntz(3)

      collect = collect - be*prj1 - prj2 - bec*prj3

c----
c factor of 1/2 from trapezoidal rule
c - sign from unfortunate definition:
c bv = n x t instead of b = t x n
c----

      A(k,l) = - 0.5000 * collect

      A(k,l) = A(k,l)/arel(k)

  20  Continue         ! run over the elements

       cel_imp(l) = 0.0D0  ! reset

  10  Continue      ! run over impulses

c----------------
c----------------
c End of Impulses
c----------------
c----------------

c--------------------------------------
c complete the finite-difference matrix
c--------------------------------------

      Dth = 0.50D0*Dt

      Do i=1,Nelm

       If(Move.eq.0) then
          A(i,i) = - Diltel(i) + A(i,i)
       Else If(Move.eq.1) then
          A(i,i) = - Dexpel(i) + A(i,i)
       End If

       Do j=1,Nelm
        A(i,j)=Dth*A(i,j)
       End Do

      End Do

c------------------------------------------
c compute the RHS of the transport equation
c------------------------------------------

      Do i=1,Nelm
       RHS(i) = 0.0D0
       Do j=1,Nelm
        RHS(i) = RHS(i) + A(i,j)*cel(j)
       End Do
       RHS(i) = cel(i) + RHS(i)
      End Do

c---------------
c compute the coefficient matrix
c of the transport equation
c---------------

      Do i=1,Nelm
       Do j=1,Nelm
        AM(i,j) = - A(i,j)
       End Do
       AM(i,i) = AM(i,i)+1.0D0
      End Do

c------------------------
c Solve the linear system
c------------------------

      Isym_gel   = 0
      Iwlpvt_gel = 1

      call gel
     +
     + (Nelm
     + ,AM
     + ,RHS
     + ,cel
     + ,Isym_gel
     + ,Iwlpvt_gel
c    + ,l,u
c    + ,det
     + ,Istop
     + )

c----------------------------------
c Display the element concentrations
c----------------------------------

      Do i=1,Nelm
c      write (6,100) i,cel(i),RHS(i)
      End Do
c     pause

c----------------------------------------------
c node concentration from element concentration
c----------------------------------------------

       call surfact_fvm_etn
     +
     +  (Npts
     +  ,cel
     +  ,c
     +  )

c-----
c done
c-----

 100  Format (1x,i5,10(1x,f15.10))
 101  Format (2(1x,i5),10(1x,f15.10))

      return
      end
