      subroutine caps_3d_slp 
     +
     +   (nelm
     +   ,npts
     +   ,mint
     +   ,NGL
     +   ,Elst
     +   ,Elstb
     +   ,crvmr
     +   ,Iflow
     +   ,slp
     +   )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------
c Computes the single-layer potential
c at the nodes of a triangular grid
c
c SYMBOLS:
c -------
c
c slp(k,i):   ith component of slp at
c             kth node
c
c srtn:   surface tension
c sgsrtn: surface gradient of surface tension
c
c crvm:  mean curvature at the nodes
c crvt:  curvature tensor at the nodes
c
c------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)
      Dimension  vna(1026,3)
      Dimension crvm(1026)
      Dimension crvt(1026,3,3)

      Dimension   nvel(1026)
      Dimension   srtn(1026)
      Dimension sgsrtn(1026,3)

      Dimension  elten(1026,3,3)     ! in plane elastic tensions
      Dimension     tst(1026,3)      ! transverse shear tension

      Dimension pxx(1026),pxy(1026),pxz(1026)
      Dimension pyx(1026),pyy(1026),pyz(1026)
      Dimension pzx(1026),pzy(1026),pzz(1026)

      Dimension Df(1026,3),slp(1026,3)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)
      Dimension  Dfel(512,4)

      Dimension  v_aux(1026,3)    ! auxiliary vector
      Dimension    vtg(1026,3,3)  ! tangential gradient of a vector

      Dimension jxy(100),arl(0:101),tnx(100),tny(100)
      Dimension eltenxy(100),eltenz(100),tstxy(100)
      Dimension angle(100)
      Dimension fsm(-10:100)

      Dimension zz(20),ww(20)
      Dimension xiq(20),etq(20),wq(20)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe

      common/albega/alpha,beta,gamma

      common/geo2/vna
      common/geo3/crvm

      common/geo8/nxy,jxy

      common/tension/srtn

      common/tenten/elten,tst    ! to be passed to Dfel

      common/var/shrt,wall

      common/veloc1/nvelt,nvel

      common/dff/Df
      common/dffel/Dfel

      common/zwl/zz,ww
      common/trq/xiq,etq,wq
      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c----------
c  Go to 99    ! to bypass the computation of the slp
c----------

c--------------
c initialize Df
c--------------

      Do i=1,npts
        Df(i,1) = 0.0
        Df(i,2) = 0.0
        Df(i,3) = 0.0
      End Do

c-----------------------------------------------
c Compute Df at the nodes due to surface tension
c-----------------------------------------------

      Go to 98    ! bypass

c-----
c compute the mean curvature at the nodes
c-----

      call crvm_3d (nelm,npts)

c---
c compute the surface gradient
c of the surface tension
c---

      call sgrad_3d
     +
     +   (npts
     +   ,nelm
     +   ,srtn
     +   ,sgsrtn
     +   )

c---
c compute Df due to surface tension
c---

c     write (6,*)
c     write (6,*)  "Traction discontinuity"
c     write (6,*)  "due to surface tension, curvature"
c     write (6,*)

      Do i=1,npts

       cf = 2.0*crvm(i)*srtn(i)

       Df(i,1) = cf*vna(i,1) - sgsrtn(i,1)
       Df(i,2) = cf*vna(i,2) - sgsrtn(i,2)
       Df(i,3) = cf*vna(i,3) - sgsrtn(i,3)

c      prj = vna(i,1)*Df(i,1) + vna(i,2)*Df(i,2)
c    +     + vna(i,3)*Df(i,3)
c
c      write (6,100) i,Df(i,1),Df(i,2),Df(i,3),crvm(i),prj

      End Do

c------------------------------------------------
c Df at the nodes due to surface tension computed
c------------------------------------------------

  98  Continue

c----------------------------------
c Df at the nodes due to elasticity
c----------------------------------

c--------------------------------------
c compute the in-plane elastic tensions
c at the nodes
c using a shell-theory constitutive
c equation
c--------------------------------------

      Do i=1,npts
        elten(i,1,1) = 0.0
        elten(i,1,2) = 0.0
        elten(i,1,3) = 0.0
        elten(i,2,1) = 0.0
        elten(i,2,2) = 0.0
        elten(i,2,3) = 0.0
        elten(i,3,1) = 0.0
        elten(i,3,2) = 0.0
        elten(i,3,3) = 0.0
      End Do

      if(Elst.ne.0) then

      call elten_3d
     +
     +   (nelm
     +   ,npts
     +   ,Elst
     +   ,elten
     +   )

      end if

c------------------------------------------------
c proceed to compute the transverse shear tension
c------------------------------------------------

      Do i=1,npts
        crvt(i,1,1) = 0.0
        crvt(i,1,2) = 0.0
        crvt(i,1,3) = 0.0
        crvt(i,2,1) = 0.0
        crvt(i,2,2) = 0.0
        crvt(i,2,3) = 0.0
        crvt(i,3,1) = 0.0
        crvt(i,3,2) = 0.0
        crvt(i,3,3) = 0.0
      End Do

      If(Elstb.eq.0) Go to 97

c-----------------------------
c compute the curvature tensor
c as the tangential gradient
c of the normal vector
c-----------------------------

      Do i=1,npts
       v_aux(i,1) = vna(i,1)     ! auxiliary vector
       v_aux(i,2) = vna(i,2)
       v_aux(i,3) = vna(i,3)
      End Do

      call vtg_3d
     +
     +  (nelm
     +  ,npts
     +  ,v_aux
     +  ,crvt
     +  )

c-------------------------------------------
c force the curvature tensor to be symmetric
c
c and project it in the tangent plane
c-------------------------------------------

      Do i=1,npts

       tmp = 0.5*(crvt(i,1,2)+crvt(i,2,1))
       crvt(i,1,2) = tmp
       crvt(i,2,1) = tmp
       tmp = 0.5*(crvt(i,1,3)+crvt(i,3,1))
       crvt(i,1,3) = tmp
       crvt(i,3,1) = tmp
       tmp = 0.5*(crvt(i,2,3)+crvt(i,3,2))
       crvt(i,2,3) = tmp
       crvt(i,3,2) = tmp

       pxx(i) = 1.0 - vna(i,1)**2          ! projection matrix
       pxy(i) =     - vna(i,1)*vna(i,2)
       pxz(i) =     - vna(i,1)*vna(i,3)
       pyy(i) = 1.0 - vna(i,2)**2
       pyz(i) =     - vna(i,2)*vna(i,3)
       pzz(i) = 1.0 - vna(i,3)**2
       pyx(i) = pxy(i)
       pzx(i) = pxz(i)
       pzy(i) = pyz(i)

       tmpxx = crvt(i,1,1)
       tmpxy = crvt(i,1,2)
       tmpxz = crvt(i,1,3)

       tmpyx = crvt(i,2,1)
       tmpyy = crvt(i,2,2)
       tmpyz = crvt(i,2,3)

       tmpzx = crvt(i,3,1)
       tmpzy = crvt(i,3,2)
       tmpzz = crvt(i,3,3)

       crvt(i,1,1) = pxx(i)*tmpxx + pxy(i)*tmpyx + pxz(i)*tmpzx
       crvt(i,1,2) = pxx(i)*tmpxy + pxy(i)*tmpyy + pxz(i)*tmpzy
       crvt(i,1,3) = pxx(i)*tmpxz + pxy(i)*tmpyz + pxz(i)*tmpzz

       crvt(i,2,1) = pyx(i)*tmpxx + pyy(i)*tmpyx + pyz(i)*tmpzx
       crvt(i,2,2) = pyx(i)*tmpxy + pyy(i)*tmpyy + pyz(i)*tmpzy
       crvt(i,2,3) = pyx(i)*tmpxz + pyy(i)*tmpyz + pyz(i)*tmpzz

       crvt(i,3,1) = pzx(i)*tmpxx + pzy(i)*tmpyx + pzz(i)*tmpzx
       crvt(i,3,2) = pzx(i)*tmpxy + pzy(i)*tmpyy + pzz(i)*tmpzy
       crvt(i,3,3) = pzx(i)*tmpxz + pzy(i)*tmpyz + pzz(i)*tmpzz

       crv_mean = 0.5*(crvt(i,1,1)+crvt(i,2,2)+crvt(i,3,3))


c      write (6,100) i,crv_mean
c      write (6,100) i,crvt(i,1,1),crvt(i,1,2),crvt(i,1,3)
c      write (6,100) i,crvt(i,2,1),crvt(i,2,2),crvt(i,2,3)
c      write (6,100) i,crvt(i,3,1),crvt(i,3,2),crvt(i,3,3)

c-----------------------
c testing
c
       Go to 888

       crvt(i,1,1) = pxx(i)
       crvt(i,1,2) = pxy(i)
       crvt(i,1,3) = pxz(i)
       crvt(i,2,1) = pyx(i)
       crvt(i,2,2) = pyy(i)
       crvt(i,2,3) = pyz(i)
       crvt(i,3,1) = pzx(i)
       crvt(i,3,2) = pzy(i)
       crvt(i,3,3) = pzz(i)

       crv_mean = 0.5*(crvt(i,1,1)+crvt(i,2,2)+crvt(i,3,3))

       write (6,100) i,crv_mean
       write (6,100) i,crvt(i,1,1),crvt(i,1,2),crvt(i,1,3)
       write (6,100) i,crvt(i,2,1),crvt(i,2,2),crvt(i,2,3)
       write (6,100) i,crvt(i,3,1),crvt(i,3,2),crvt(i,3,3)

 888  Continue

c end of testing
c----------------

c----
c subtract the resting shape
c mean curvature
c-----

       crvt(i,1,1) = crvt(i,1,1) - crvmr * pxx(i)
       crvt(i,1,2) = crvt(i,1,2) - crvmr * pxy(i)
       crvt(i,1,3) = crvt(i,1,3) - crvmr * pxz(i)
       crvt(i,2,1) = crvt(i,2,1) - crvmr * pyx(i)
       crvt(i,2,2) = crvt(i,2,2) - crvmr * pyy(i)
       crvt(i,2,3) = crvt(i,2,3) - crvmr * pyz(i)
       crvt(i,3,1) = crvt(i,3,1) - crvmr * pzx(i)
       crvt(i,3,2) = crvt(i,3,2) - crvmr * pzy(i)
       crvt(i,3,3) = crvt(i,3,3) - crvmr * pzz(i)

      End Do

c-------------------------------------
c compute the transverse shear tension
c as the surface divergence of the
c curvature tensor
c-------------------------------------

      Do i=1,npts
       v_aux(i,1) = crvt(i,1,1)     ! auxilliary vector
       v_aux(i,2) = crvt(i,2,1)
       v_aux(i,3) = crvt(i,3,1)
      End Do

      call vtg_3d
     +
     +  (nelm
     +  ,npts
     +  ,v_aux
     +  ,vtg
     +  )

c---
c pick up the trace
c---

      Do i=1,npts
       tst(i,1) = vtg(i,1,1)+vtg(i,2,2)+vtg(i,3,3)
      End Do

      Do i=1,npts
       v_aux(i,1) = crvt(i,1,2)     ! auxilliary vector
       v_aux(i,2) = crvt(i,2,2)
       v_aux(i,3) = crvt(i,3,2)
      End Do

      call vtg_3d
     +
     +   (nelm
     +   ,npts
     +   ,v_aux
     +   ,vtg
     +   )

c---
c pick up the trace
c---

      Do i=1,npts
       tst(i,2) = vtg(i,1,1)+vtg(i,2,2)+vtg(i,3,3)
      End Do

      Do i=1,npts
       v_aux(i,1) = crvt(i,1,3)     ! auxilliary vector
       v_aux(i,2) = crvt(i,2,3)
       v_aux(i,3) = crvt(i,3,3)
      End Do

      call vtg_3d
     +
     +   (nelm
     +   ,npts
     +   ,v_aux
     +   ,vtg
     +   )

c---
c pick up the trace
c---

      Do i=1,npts
       tst(i,3) = vtg(i,1,1)+vtg(i,2,2)+vtg(i,3,3)
      End Do

c------------------------
c project tst onto P=I-nn
c------------------------

c     write (6,*)
c     write (6,*) "Mean curvature, transverse shear tensions"
c     write (6,*)

      tstmax = 0.0

      Do i=1,npts

       tmpx = tst(i,1)
       tmpy = tst(i,2)
       tmpz = tst(i,3)

       tst(i,1) = Elstb*(pxx(i)*tmpx + pxy(i)*tmpy + pxz(i)*tmpz)
       tst(i,2) = Elstb*(pyx(i)*tmpx + pyy(i)*tmpy + pyz(i)*tmpz)
       tst(i,3) = Elstb*(pzx(i)*tmpx + pzy(i)*tmpy + pzz(i)*tmpz)

c--------------- testing
       Go to 887

       If(ne(i,1).ne.2) then
       tstm = sqrt(tst(i,1)**2+tst(i,2)**2+tst(i,3)**2)
       tstn = tst(i,1)*vna(i,1)+tst(i,2)*vna(i,2)
     +           +tst(i,3)*vna(i,3)
       If(tstm.gt.tstmax) tstmax = tstm
       crv_mean = 0.5*(crvt(i,1,1)+crvt(i,2,2)+crvt(i,3,3))

       write (6,100) i,crv_mean,tst(i,1),tst(i,2),tst(i,3),tstm
     +               ,tstn
       End If

 887   Continue
c--------------

      End Do

c     write (6,109) tstmax

  97  Continue

c-----------------------
c printing elten and tst
c-----------------------

      Go to 886

      open (7,file="PLOTDAT")

      write (7,100) nxy

c---
c compute the polygonal arc length
c---

      arl(1) = 0.0

      Do i=2,nxy
       ia = i-1
       m  = jxy(i)
       ma = jxy(ia)
       arl(i) = arl(ia) + sqrt((p(m,1)-p(ma,1))**2
     +                        +(p(m,2)-p(ma,2))**2)
      End Do

      arl(0) = -arl(2)
      arl(nxy+1) = arl(nxy)+arl(2)

c----
c compute the tensions
c---

      Do i=1,nxy

       ia = i-1
       i1 = i+1

       ma = jxy(ia)
       m  = jxy(i)
       m1 = jxy(i1)

c---
c compute the unit tangent vector
c---

       x0 = arl(ia)-arl(i)
       x1 = arl(i1)-arl(i)

       If(i.eq.1) then
        y0 = p(jxy(nxy-1),1)-p(m,1)
       Else
        y0 = p(ma,1)-p(m,1)
       End If
       y1 = p(m1,1)-p(m,1)

       DxDl = (x0*y1/x1 - x1*y0/x0)/(x0-x1)

       If(i.eq.1) then
        y0 = p(jxy(nxy-1),2)-p(m,2)
       Else
        y0 = p(ma,2)-p(m,2)
       End If
       y1 = p(m1,2)-p(m,2)

       DyDl = (x0*y1/x1 - x1*y0/x0)/(x0-x1)

       DlDl = sqrt(DxDl**2+DyDl**2)   ! normalize
       tnx(i) = DxDl/DlDl
       tny(i) = DyDl/DlDl

       If(i.eq.nxy) then
        tnx(i) = tnx(1)
        tny(i) = tny(1)
       End If

c---
c compute the inclination angle
c---

       If(i.eq.1) then
          xxx      = p(m,1)
          yyy      = p(m,2)
          angle(1) = atan2(yyy,xxx)
       Else
          xxx1 = p(m,1)
          yyy1 = p(m,2)
          xxx2 = p(ma,1)
          yyy2 = p(ma,2)
          prj  = xxx1*xxx2+yyy1*yyy2
          rr1  = sqrt(xxx1**2+yyy1**2)
          rr2  = sqrt(xxx2**2+yyy2**2)
          arg  = prj/(rr1*rr2)
          dan  = acos(arg)
          angle(i) = angle(ia)-dan
       End If

c      tstxy(i) = sqrt(tst(m,1)**2+tst(m,2)**2+tst(m,3)**2)
c      eltenxy(i) = sqrt(elten(m,1,1)**2+elten(m,1,2)**2)

       tstxy(i) = tst(m,1)*tnx(i)+tst(m,2)*tny(i)
       eltenx = elten(m,1,1)*tnx(i)+elten(m,1,2)*tny(i)
       elteny = elten(m,2,1)*tnx(i)+elten(m,2,2)*tny(i)
       eltenxy(i) = eltenx*tnx(i)+elteny*tny(i)
       eltenz (i) = elten(m,3,3)

      End Do

c---
c smooth tstxy
c---

      Nsm = 5

      Do Ism=1,Nsm

      Do i=1,nxy
       fsm(i) = tstxy(i)
      End Do
      fsm( 0) = fsm(nxy-1)
      fsm(-1) = fsm(nxy-2)
      fsm(nxy+1) = fsm(2)
      fsm(nxy+2) = fsm(3)

      Do i=1,nxy
       tstxy(i) = (-fsm(i-2)+4.0*fsm(i-1)+10.0*fsm(i)
     +                      +4.0*fsm(i+1)-fsm(i+2))/16.0
      End Do

c---
c smoooth eltenxy
c---

      Do i=1,nxy
       fsm(i) = eltenxy(i)
      End Do

      fsm( 0) = fsm(nxy-1)
      fsm(-1) = fsm(nxy-2)
      fsm(nxy+1) = fsm(2)
      fsm(nxy+2) = fsm(3)

      Do i=1,nxy
       eltenxy(i) = (-fsm(i-2)+4.0*fsm(i-1)+10.0*fsm(i)
     +                        +4.0*fsm(i+1)-fsm(i+2))/16.0
      End Do

c---
c smoooth eltenz
c---

      Do i=1,nxy
       fsm(i) = eltenz(i)
      End Do
      fsm( 0) = fsm(nxy-1)
      fsm(-1) = fsm(nxy-2)
      fsm(nxy+1) = fsm(2)
      fsm(nxy+2) = fsm(3)

      Do i=1,nxy
       eltenz(i) = (-fsm(i-2)+4.0*fsm(i-1)+10.0*fsm(i)
     +                       +4.0*fsm(i+1)-fsm(i+2))/16.0
      End Do

      End Do

c------
c printing
c------

      Do i=1,nxy
       prangle = 2.0+angle(i)/pi
       m  = jxy(i)
c        write (6,100) i,prangle,p(m,1),p(m,2),p(m,3)
c        write (6,101) elten(m,1,1),elten(m,1,2),elten(m,1,3)
c        write (6,101) elten(m,2,1),elten(m,2,2),elten(m,2,3)
c        write (6,101) elten(m,3,1),elten(m,3,2),elten(m,3,3)
c        write (6,101) tst(m,1),tst(m,2),tst(m,3)
         write (6,100) i,prangle,eltenxy(i),eltenz(i),tstxy(i)
         write (7,100) i,prangle,eltenxy(i),eltenz(i),tstxy(i)
      End Do

      null = 0
      write (7,100) null
      close (7)
      stop

 886  Continue

c----------------------------------------
c compute elastic traction discontinuity
c integrated over the elements
c----------------------------------------

      call caps_3d_dfel
     +
     +   (npts
     +   ,nelm
     +   )

c---------
c printing
c---------
c
c     write (6,*)
c     write (6,*) " Df averaged over the elements:"
c     write (6,*)
c
c     Do k=1,nelm
c      write (6,100) k,(Dfel(k,l),l=1,4)
c     End Do

c---------------------------
c Proceed to compute the slp
c---------------------------

c-------------------------------------
c Loop over velocity evaluation points
c-------------------------------------

      Do 1 node=1,nvelt

       i = nvel(node)  ! i is the global node label

c      write (6,*) " Computing the slp at point",i

       x0 = p(i,1)
       y0 = p(i,2)
       z0 = p(i,3)

       us = 0.0
       vs = 0.0
       ws = 0.0

       Do 2 k=1,nelm     ! run over elements

        uxel = 0.0
        uyel = 0.0
        uzel = 0.0

        if(    i.ne.n(k,1).and.i.ne.n(k,2)
     +    .and.i.ne.n(k,3).and.i.ne.n(k,4)
     +    .and.i.ne.n(k,5).and.i.ne.n(k,6)
     +    ) then

c----------------------------------
c Non-singular integration:
c
c use a regular triangle quadrature
c----------------------------------

        call caps_3d_slp_integrate
     +
     +    (x0,y0,z0
     +    ,k           ! element index
     +    ,mint
     +    ,Iflow
     +    ,uxel,uyel,uzel
     +    )

        Go to 3                ! Do another element

        End If

c-------------------------------------------
c Singular integration:
c
c If the point i is a vertex node,
c will integrate over the flat triangle
c defined by the node
c using the polar integration rule
c
c If the point i is a mid node,
c will breakup the curved triangle into four
c flat triangles
c and integrate over the flat triangles
c using the polar integration rule
c  
c   Iopt_int = 1 only the position vector
c              2 position vector and rest of variables
c
c--------------------------------------------

      i1 = n(k,1)  ! global node index

      i2 = n(k,2)
      i3 = n(k,3)
      i4 = n(k,4)
      i5 = n(k,5)
      i6 = n(k,6)

c--------------------------------------------
c singular element with singularity at node 1
c vertex node
c--------------------------------------------

       If(i.eq.n(k,1)) then

          x1 =  p(i1,1)
          y1 =  p(i1,2)
          z1 =  p(i1,3)

         fx1 = Df(i1,1)+Dfel(k,1)
         fy1 = Df(i1,2)+Dfel(k,2)
         fz1 = Df(i1,3)+Dfel(k,3)

          x2 =  p(i2,1)
          y2 =  p(i2,2)
          z2 =  p(i2,3)

         fx2 = Df(i2,1)+Dfel(k,1)
         fy2 = Df(i2,2)+Dfel(k,2)
         fz2 = Df(i2,3)+Dfel(k,3)

          x3 =  p(i3,1)
          y3 =  p(i3,2)
          z3 =  p(i3,3)

         fx3 = Df(i3,1)+Dfel(k,1)
         fy3 = Df(i3,2)+Dfel(k,2)
         fz3 = Df(i3,3)+Dfel(k,3)

         call caps_3d_slp_sing
     +
     +      (NGL,Iflow
     +      ,x1,y1,z1
     +      ,x2,y2,z2
     +      ,x3,y3,z3
     +      ,fx1,fy1,fz1
     +      ,fx2,fy2,fz2
     +      ,fx3,fy3,fz3
     +      ,uxel,uyel,uzel
     +      ,GExx,GExy,GExz
     +      ,GEyx,GEyy,GEyz
     +      ,GEzx,GEzy,GEzz 
     +      )

c--------------------------------------------
c singular element with singularity at node 2
c vertex node
c--------------------------------------------

        Else If(i.eq.n(k,2)) then

          x1 =  p(i2,1)
          y1 =  p(i2,2)
          z1 =  p(i2,3)
         fx1 = Df(i2,1)+Dfel(k,1)
         fy1 = Df(i2,2)+Dfel(k,2)
         fz1 = Df(i2,3)+Dfel(k,3)

          x2 =  p(i3,1)
          y2 =  p(i3,2)
          z2 =  p(i3,3)
         fx2 = Df(i3,1)+Dfel(k,1)
         fy2 = Df(i3,2)+Dfel(k,2)
         fz2 = Df(i3,3)+Dfel(k,3)

          x3 =  p(i1,1)
          y3 =  p(i1,2)
          z3 =  p(i1,3)
         fx3 = Df(i1,1)+Dfel(k,1)
         fy3 = Df(i1,2)+Dfel(k,2)
         fz3 = Df(i1,3)+Dfel(k,3)

         call caps_3d_slp_sing
     +
     +       (NGL,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,fx1,fy1,fz1
     +       ,fx2,fy2,fz2
     +       ,fx3,fy3,fz3
     +       ,uxel,uyel,uzel
     +       ,GExx,GExy,GExz
     +       ,GEyx,GEyy,GEyz
     +       ,GEzx,GEzy,GEzz 
     +       )

c--------------------------------------------
c singular element with singularity at node 3
c vertex node
c--------------------------------------------

        else if(i.eq.n(k,3)) then

          x1 =  p(i3,1)
          y1 =  p(i3,2)
          z1 =  p(i3,3)
         fx1 = Df(i3,1)+Dfel(k,1)
         fy1 = Df(i3,2)+Dfel(k,2)
         fz1 = Df(i3,3)+Dfel(k,3)

          x2 =  p(i1,1)
          y2 =  p(i1,2)
          z2 =  p(i1,3)
         fx2 = Df(i1,1)+Dfel(k,1)
         fy2 = Df(i1,2)+Dfel(k,2)
         fz2 = Df(i1,3)+Dfel(k,3)

          x3 =  p(i2,1)
          y3 =  p(i2,2)
          z3 =  p(i2,3)
         fx3 = Df(i2,1)+Dfel(k,1)
         fy3 = Df(i2,2)+Dfel(k,2)
         fz3 = Df(i2,3)+Dfel(k,3)

         call caps_3d_slp_sing
     +
     +       (NGL,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,fx1,fy1,fz1
     +       ,fx2,fy2,fz2
     +       ,fx3,fy3,fz3
     +       ,uxel,uyel,uzel
     +       ,GExx,GExy,GExz
     +       ,GEyx,GEyy,GEyz
     +       ,GEzx,GEzy,GEzz 
     +       )

c--------------------------------------------
c singular element with singularity at node 4
c mid-node
c Will integrate over 4 flat triangles
c--------------------------------------------

        else if(i.eq.n(k,4)) then

          x1 =  p(i4,1)
          y1 =  p(i4,2)
          z1 =  p(i4,3)
         fx1 = Df(i4,1)+Dfel(k,1)
         fy1 = Df(i4,2)+Dfel(k,2)
         fz1 = Df(i4,3)+Dfel(k,3)

          x2 =  p(i6,1)
          y2 =  p(i6,2)
          z2 =  p(i6,3)
         fx2 = Df(i6,1)+Dfel(k,1)
         fy2 = Df(i6,2)+Dfel(k,2)
         fz2 = Df(i6,3)+Dfel(k,3)

          x3 =  p(i1,1)
          y3 =  p(i1,2)
          z3 =  p(i1,3)
         fx3 = Df(i1,1)+Dfel(k,1)
         fy3 = Df(i1,2)+Dfel(k,2)
         fz3 = Df(i1,3)+Dfel(k,3)

         call caps_3d_slp_sing
     +
     +       (NGL,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,fx1,fy1,fz1
     +       ,fx2,fy2,fz2
     +       ,fx3,fy3,fz3
     +       ,uxel,uyel,uzel
     +       ,GExx,GExy,GExz
     +       ,GEyx,GEyy,GEyz
     +       ,GEzx,GEzy,GEzz 
     +       )

          x1 =  p(i4,1)
          y1 =  p(i4,2)
          z1 =  p(i4,3)
         fx1 = Df(i4,1)+Dfel(k,1)
         fy1 = Df(i4,2)+Dfel(k,2)
         fz1 = Df(i4,3)+Dfel(k,3)

          x2 =  p(i3,1)
          y2 =  p(i3,2)
          z2 =  p(i3,3)
         fx2 = Df(i3,1)+Dfel(k,1)
         fy2 = Df(i3,2)+Dfel(k,2)
         fz2 = Df(i3,3)+Dfel(k,3)

          x3 =  p(i6,1)
          y3 =  p(i6,2)
          z3 =  p(i6,3)
         fx3 = Df(i6,1)+Dfel(k,1)
         fy3 = Df(i6,2)+Dfel(k,2)
         fz3 = Df(i6,3)+Dfel(k,3)

         call caps_3d_slp_sing
     +
     +       (NGL,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,fx1,fy1,fz1
     +       ,fx2,fy2,fz2
     +       ,fx3,fy3,fz3
     +       ,uxel,uyel,uzel
     +       ,GExx,GExy,GExz
     +       ,GEyx,GEyy,GEyz
     +       ,GEzx,GEzy,GEzz 
     +       )

          x1 =  p(i4,1)
          y1 =  p(i4,2)
          z1 =  p(i4,3)
         fx1 = Df(i4,1)+Dfel(k,1)
         fy1 = Df(i4,2)+Dfel(k,2)
         fz1 = Df(i4,3)+Dfel(k,3)

          x2 =  p(i5,1)
          y2 =  p(i5,2)
          z2 =  p(i5,3)
         fx2 = Df(i5,1)+Dfel(k,1)
         fy2 = Df(i5,2)+Dfel(k,2)
         fz2 = Df(i5,3)+Dfel(k,3)

          x3 =  p(i3,1)
          y3 =  p(i3,2)
          z3 =  p(i3,3)
         fx3 = Df(i3,1)+Dfel(k,1)
         fy3 = Df(i3,2)+Dfel(k,2)
         fz3 = Df(i3,3)+Dfel(k,3)

         call caps_3d_slp_sing
     +
     +         (NGL,Iflow
     +         ,x1,y1,z1
     +         ,x2,y2,z2
     +         ,x3,y3,z3
     +         ,fx1,fy1,fz1
     +         ,fx2,fy2,fz2
     +         ,fx3,fy3,fz3
     +         ,uxel,uyel,uzel
     +         ,GExx,GExy,GExz
     +         ,GEyx,GEyy,GEyz
     +         ,GEzx,GEzy,GEzz
     +         )

          x1 =  p(i4,1)
          y1 =  p(i4,2)
          z1 =  p(i4,3)
         fx1 = Df(i4,1)+Dfel(k,1)
         fy1 = Df(i4,2)+Dfel(k,2)
         fz1 = Df(i4,3)+Dfel(k,3)

          x2 =  p(i2,1)
          y2 =  p(i2,2)
          z2 =  p(i2,3)
         fx2 = Df(i2,1)+Dfel(k,1)
         fy2 = Df(i2,2)+Dfel(k,2)
         fz2 = Df(i2,3)+Dfel(k,3)

          x3 =  p(i5,1)
          y3 =  p(i5,2)
          z3 =  p(i5,3)
         fx3 = Df(i5,1)+Dfel(k,1)
         fy3 = Df(i5,2)+Dfel(k,2)
         fz3 = Df(i5,3)+Dfel(k,3)

         call caps_3d_slp_sing
     +
     +        (NGL,Iflow
     +        ,x1,y1,z1
     +        ,x2,y2,z2
     +        ,x3,y3,z3
     +        ,fx1,fy1,fz1
     +        ,fx2,fy2,fz2
     +        ,fx3,fy3,fz3
     +        ,uxel,uyel,uzel
     +        ,GExx,GExy,GExz
     +        ,GEyx,GEyy,GEyz
     +        ,GEzx,GEzy,GEzz
     +        )

c--------------------------------------------
c singular element with singularity at node 5
c mid-node
c Will integrate over 4 flat triangles
c--------------------------------------------

        Else If(i.eq.n(k,5)) then

          x1 =  p(i5,1)
          y1 =  p(i5,2)
          z1 =  p(i5,3)
         fx1 = Df(i5,1)+Dfel(k,1)
         fy1 = Df(i5,2)+Dfel(k,2)
         fz1 = Df(i5,3)+Dfel(k,3)

          x2 =  p(i4,1)
          y2 =  p(i4,2)
          z2 =  p(i4,3)
         fx2 = Df(i4,1)+Dfel(k,1)
         fy2 = Df(i4,2)+Dfel(k,2)
         fz2 = Df(i4,3)+Dfel(k,3)

          x3 =  p(i2,1)
          y3 =  p(i2,2)
          z3 =  p(i2,3)
         fx3 = Df(i2,1)+Dfel(k,1)
         fy3 = Df(i2,2)+Dfel(k,2)
         fz3 = Df(i2,3)+Dfel(k,3)

         call caps_3d_slp_sing
     +
     +       (NGL,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,fx1,fy1,fz1
     +       ,fx2,fy2,fz2
     +       ,fx3,fy3,fz3
     +       ,uxel,uyel,uzel
     +       ,GExx,GExy,GExz
     +       ,GEyx,GEyy,GEyz
     +       ,GEzx,GEzy,GEzz
     +       )

          x1 =  p(i5,1)
          y1 =  p(i5,2)
          z1 =  p(i5,3)
         fx1 = Df(i5,1)+Dfel(k,1)
         fy1 = Df(i5,2)+Dfel(k,2)
         fz1 = Df(i5,3)+Dfel(k,3)

          x2 =  p(i1,1)
          y2 =  p(i1,2)
          z2 =  p(i1,3)
         fx2 = Df(i1,1)+Dfel(k,1)
         fy2 = Df(i1,2)+Dfel(k,2)
         fz2 = Df(i1,3)+Dfel(k,3)

          x3 =  p(i4,1)
          y3 =  p(i4,2)
          z3 =  p(i4,3)
         fx3 = Df(i4,1)+Dfel(k,1)
         fy3 = Df(i4,2)+Dfel(k,2)
         fz3 = Df(i4,3)+Dfel(k,3)

         call caps_3d_slp_sing
     +
     +        (NGL,Iflow
     +        ,x1,y1,z1
     +        ,x2,y2,z2
     +        ,x3,y3,z3
     +        ,fx1,fy1,fz1
     +        ,fx2,fy2,fz2
     +        ,fx3,fy3,fz3
     +        ,uxel,uyel,uzel
     +        ,GExx,GExy,GExz
     +        ,GEyx,GEyy,GEyz
     +        ,GEzx,GEzy,GEzz
     +        )

          x1 =  p(i5,1)
          y1 =  p(i5,2)
          z1 =  p(i5,3)
         fx1 = Df(i5,1)+Dfel(k,1)
         fy1 = Df(i5,2)+Dfel(k,2)
         fz1 = Df(i5,3)+Dfel(k,3)

          x2 =  p(i6,1)
          y2 =  p(i6,2)
          z2 =  p(i6,3)
         fx2 = Df(i6,1)+Dfel(k,1)
         fy2 = Df(i6,2)+Dfel(k,2)
         fz2 = Df(i6,3)+Dfel(k,3)

          x3 =  p(i1,1)
          y3 =  p(i1,2)
          z3 =  p(i1,3)
         fx3 = Df(i1,1)+Dfel(k,1)
         fy3 = Df(i1,2)+Dfel(k,2)
         fz3 = Df(i1,3)+Dfel(k,3)

         call caps_3d_slp_sing
     +
     +        (NGL,Iflow
     +        ,x1,y1,z1
     +        ,x2,y2,z2
     +        ,x3,y3,z3
     +        ,fx1,fy1,fz1
     +        ,fx2,fy2,fz2
     +        ,fx3,fy3,fz3
     +        ,uxel,uyel,uzel
     +        ,GExx,GExy,GExz
     +        ,GEyx,GEyy,GEyz
     +        ,GEzx,GEzy,GEzz
     +        )

          x1 =  p(i5,1)
          y1 =  p(i5,2)
          z1 =  p(i5,3)
         fx1 = Df(i5,1)+Dfel(k,1)
         fy1 = Df(i5,2)+Dfel(k,2)
         fz1 = Df(i5,3)+Dfel(k,3)

          x2 =  p(i3,1)
          y2 =  p(i3,2)
          z2 =  p(i3,3)
         fx2 = Df(i3,1)+Dfel(k,1)
         fy2 = Df(i3,2)+Dfel(k,2)
         fz2 = Df(i3,3)+Dfel(k,3)

          x3 =  p(i6,1)
          y3 =  p(i6,2)
          z3 =  p(i6,3)
         fx3 = Df(i6,1)+Dfel(k,1)
         fy3 = Df(i6,2)+Dfel(k,2)
         fz3 = Df(i6,3)+Dfel(k,3)

         call caps_3d_slp_sing
     +
     +         (NGL,Iflow
     +         ,x1,y1,z1
     +         ,x2,y2,z2
     +         ,x3,y3,z3
     +         ,fx1,fy1,fz1
     +         ,fx2,fy2,fz2
     +         ,fx3,fy3,fz3
     +         ,uxel,uyel,uzel
     +         ,GExx,GExy,GExz
     +         ,GEyx,GEyy,GEyz
     +         ,GEzx,GEzy,GEzz
     +         )

c--------------------------------------------
c singular element with singularity at node 6
c mid-node
c Will integrate over 4 flat triangles
c--------------------------------------------

        Else If(i.eq.n(k,6)) then

          x1 =  p(i6,1)
          y1 =  p(i6,2)
          z1 =  p(i6,3)
         fx1 = Df(i6,1)+Dfel(k,1)
         fy1 = Df(i6,2)+Dfel(k,2)
         fz1 = Df(i6,3)+Dfel(k,3)

          x2 =  p(i1,1)
          y2 =  p(i1,2)
          z2 =  p(i1,3)
         fx2 = Df(i1,1)+Dfel(k,1)
         fy2 = Df(i1,2)+Dfel(k,2)
         fz2 = Df(i1,3)+Dfel(k,3)

          x3 =  p(i4,1)
          y3 =  p(i4,2)
          z3 =  p(i4,3)
         fx3 = Df(i4,1)+Dfel(k,1)
         fy3 = Df(i4,2)+Dfel(k,2)
         fz3 = Df(i4,3)+Dfel(k,3)

         call caps_3d_slp_sing
     +
     +        (NGL,Iflow
     +        ,x1,y1,z1
     +        ,x2,y2,z2
     +        ,x3,y3,z3
     +        ,fx1,fy1,fz1
     +        ,fx2,fy2,fz2
     +        ,fx3,fy3,fz3
     +        ,uxel,uyel,uzel
     +        ,GExx,GExy,GExz
     +        ,GEyx,GEyy,GEyz
     +        ,GEzx,GEzy,GEzz
     +        )

          x1 =  p(i6,1)
          y1 =  p(i6,2)
          z1 =  p(i6,3)
         fx1 = Df(i6,1)+Dfel(k,1)
         fy1 = Df(i6,2)+Dfel(k,2)
         fz1 = Df(i6,3)+Dfel(k,3)

          x2 =  p(i4,1)
          y2 =  p(i4,2)
          z2 =  p(i4,3)
         fx2 = Df(i4,1)+Dfel(k,1)
         fy2 = Df(i4,2)+Dfel(k,2)
         fz2 = Df(i4,3)+Dfel(k,3)

          x3 =  p(i2,1)
          y3 =  p(i2,2)
          z3 =  p(i2,3)
         fx3 = Df(i2,1)+Dfel(k,1)
         fy3 = Df(i2,2)+Dfel(k,2)
         fz3 = Df(i2,3)+Dfel(k,3)

         call caps_3d_slp_sing
     +
     +       (NGL,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,fx1,fy1,fz1
     +       ,fx2,fy2,fz2
     +       ,fx3,fy3,fz3
     +       ,uxel,uyel,uzel
     +       ,GExx,GExy,GExz
     +       ,GEyx,GEyy,GEyz
     +       ,GEzx,GEzy,GEzz
     +       )

          x1 =  p(i6,1)
          y1 =  p(i6,2)
          z1 =  p(i6,3)
         fx1 = Df(i6,1)+Dfel(k,1)
         fy1 = Df(i6,2)+Dfel(k,2)
         fz1 = Df(i6,3)+Dfel(k,3)

          x2 =  p(i2,1)
          y2 =  p(i2,2)
          z2 =  p(i2,3)
         fx2 = Df(i2,1)+Dfel(k,1)
         fy2 = Df(i2,2)+Dfel(k,2)
         fz2 = Df(i2,3)+Dfel(k,3)

          x3 =  p(i5,1)
          y3 =  p(i5,2)
          z3 =  p(i5,3)
         fx3 = Df(i5,1)+Dfel(k,1)
         fy3 = Df(i5,2)+Dfel(k,2)
         fz3 = Df(i5,3)+Dfel(k,3)

         call caps_3d_slp_sing
     +
     +       (NGL,Iflow
     +       ,x1,y1,z1
     +       ,x2,y2,z2
     +       ,x3,y3,z3
     +       ,fx1,fy1,fz1
     +       ,fx2,fy2,fz2
     +       ,fx3,fy3,fz3
     +       ,uxel,uyel,uzel
     +       ,GExx,GExy,GExz
     +       ,GEyx,GEyy,GEyz
     +       ,GEzx,GEzy,GEzz
     +       )

          x1 =  p(i6,1)
          y1 =  p(i6,2)
          z1 =  p(i6,3)
         fx1 = Df(i6,1)+Dfel(k,1)
         fy1 = Df(i6,2)+Dfel(k,2)
         fz1 = Df(i6,3)+Dfel(k,3)

          x2 =  p(i5,1)
          y2 =  p(i5,2)
          z2 =  p(i5,3)
         fx2 = Df(i5,1)+Dfel(k,1)
         fy2 = Df(i5,2)+Dfel(k,2)
         fz2 = Df(i5,3)+Dfel(k,3)

          x3 =  p(i3,1)
          y3 =  p(i3,2)
          z3 =  p(i3,3)
         fx3 = Df(i3,1)+Dfel(k,1)
         fy3 = Df(i3,2)+Dfel(k,2)
         fz3 = Df(i3,3)+Dfel(k,3)

         call caps_3d_slp_sing
     +
     +         (NGL,Iflow
     +         ,x1,y1,z1
     +         ,x2,y2,z2
     +         ,x3,y3,z3
     +         ,fx1,fy1,fz1
     +         ,fx2,fy2,fz2
     +         ,fx3,fy3,fz3
     +         ,uxel,uyel,uzel
     +         ,GExx,GExy,GExz
     +         ,GEyx,GEyy,GEyz
     +         ,GEzx,GEzy,GEzz
     +         )

c------------
       end if    ! done integrating over a singular element
c------------

  3    Continue

       us = us + uxel
       vs = vs + uyel
       ws = ws + uzel

  2   Continue

      slp(i,1) = - us/pi8
      slp(i,2) = - vs/pi8
      slp(i,3) = - ws/pi8

  1   Continue               ! loop over nodes

c-----
c Done
c-----

  99  Continue

 100  Format (1x,i4,10(1x,f12.8))
 101  Format (10(1x,f12.8))
 109  Format (f15.10)

      Return
      End
