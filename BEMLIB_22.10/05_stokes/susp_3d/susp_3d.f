      program susp_3d 

c=========================================
c FEMLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c------------------------------------------
c Dynamic simulation of a suspension of
c freely-suspended particles
c for several flow configurations
c
c SYMBOLS:
c -------
c
c Nprtcl: number of particles
c
c p: particle geometrical nodes
c n: geometrical connectivity matrix
c
c pg(k,3,l): global-global geometrical nodes
c            of the lth particle
c
c pcl(i,3):  collocation nodes
c dpl(i,3):  dipole
c
c uclp(i,3):  mutually induced particle velocity
c
c Implemented capacity: 2 particles
c---------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension p(1026,3),pg(1026,3,2)
      Dimension n(512,6)

      Dimension       pcl(2306,3),pclg(2306,3,2) 
      Dimension       pev(2306,3)
      Dimension       dpl(2306,3),dplg(2306,3,2)
      Dimension   dpl_bef(2306,3)
      Dimension      dplc(2306,3)
      Dimension      uclp(2306,3),uclpg(2306,3,2),uclpp(2306,3,2,2)

      Dimension stress(3,3),stresslet(3,3)

      Dimension xiq(20),etq(20),wq(20)

      Dimension dpl1(2306,3)
      Dimension dpl2(2306,3)
      Dimension dpl3(2306,3)
      Dimension dpl4(2306,3)

c----------
c particles (up to 2)
c----------

      Dimension boa(2),coa(2),req(2)
      Dimension cx(2),cy(2),cz(2)
      Dimension dpldn(2)
      Dimension Vx(2),Vy(2),Vz(2)
      Dimension Ox(2),Oy(2),Oz(2)
      Dimension theta(2),dtheta(2),dthjef(2)
      Dimension phi(2),dphi(2)
      Dimension chi(2)
      Dimension dirx(2),diry(2),dirz(2)

      Dimension xfadh(2),yfadh(2),zfadh(2)

      Dimension Vx_old(2),Vy_old(2),Vz_old(2)
      Dimension Ox_old(2),Oy_old(2),Oz_old(2)

      Dimension cx_save(2),cy_save(2),cz_save(2)
      Dimension Vx_save(2),Vy_save(2),Vz_save(2)
      Dimension Ox_save(2),Oy_save(2),Oz_save(2)
      Dimension theta_save(2),dtheta_save(2)
      Dimension phi_save(2),dphi_save(2)

c--------------
c time stepping
c--------------

      Dimension time(5000)

c--------------
c common blocks
c--------------

      common/var1/Iflow,Iwall
      common/var2/Uinf,shrt
      common/var3/wall
      common/var4/visc

      common/var9/dpl1,dpl2,dpl3,dpl4

      common/trq/xiq,etq,wq

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c: periodic flow

      common/aaaa/a11,a12,a13,a21,a22,a23,a31,a32,a33
      common/bbbb/b11,b12,b13,b21,b22,b23,b31,b32,b33
      common/ewew/ew,tau
      common/mmmm/Max1,Max2

c----------
c constants
c----------

      pi = 3.14159 265358 D0

      piq = 0.25D0 *pi
      pih = 0.50D0 *pi
      pi2 = 2.00D0 *pi
      pi4 = 4.00D0 *pi
      pi6 = 6.00D0 *pi
      pi8 = 8.00D0 *pi
      spi4 = Dsqrt(pi4)

      null  = 0
      Nfour  = 4
      Nseven = 7
      oot  = 1.0D0/3.0D0

c-----------
c input data
c-----------

      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) '1 to read data from file: susp_3d.dat'
      write (6,*) '0 to quit'
      write (6,*) ' ---------'
      read  (5,*)  Ienrd

      if(Ienrd.eq.0) Go to 99

c read various parameters

      open (2,file="susp_3d.dat")

       read (2,*) Ioctaicos
       read (2,*) ndiv
       read (2,*)
       read (2,*) Isize
       read (2,*) intm
       read (2,*) mint
       read (2,*) tolin
       read (2,*) tolout
       read (2,*) Nter
       read (2,*) 
       read (2,*) visc
       read (2,*) 
       read (2,*) wall
       read (2,*) 
       read (2,*) Iadhere
       read (2,*) torque
       read (2,*) adhesion
       read (2,*) 
       read (2,*) Uinf
       read (2,*) shrt
       read (2,*) 
       read (2,*) Max1,Max2
       read (2,*) 
       read (2,*) mpoly
       read (2,*) Ipd      ! 0 for uniform, 1 for spectral
       read (2,*) 
       read (2,*) IRK
       read (2,*) Dt
       read (2,*) 

       if((ndiv.lt.0).or.(lt.gt.3))  then
         write (6,*) 'out of range; please try again'
         Go to 99
       end if

c read the particle position and orientation

      open (3,file="susp_3d_p.dat")

       read (3,*) Nprtcl

       Do i=1,Nprtcl
        read (3,*) 
        read (3,*) boa(i),coa(i)
        read (3,*) req(i)
        read (3,*) cx(i)
        read (3,*) cy(i)
        read (3,*) cz(i)
        read (3,*) theta(i)
        read (3,*) phi(i)
       End Do

      close(3)

c-------------------------------------
c specify the incident velocity: u_inf
c-------------------------------------

      write (6,*)
      write (6,*) 'CHOOSE THE FLOW'
      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 1 for u_x = U'
      write (6,*) ' 2 for u_y = U'
      write (6,*) ' 3 for u_z = U'
      write (6,*) ' 4 for u_x = s y'
      write (6,*) ' 5 for u_z = s y'
      write (6,*) ' 6 for u_y = s x'
      write (6,*) '10 for u_x = s (y-w)'
      write (6,*) '20 for u_x = s y, periodic flow'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) Iflow

      if(Iflow.eq.0) Go to 99

c------------------
c open output files
c------------------

      open (1,file="susp_3d.rhe")
      open (3,file="susp_3d.trj")
      open (4,file="susp_3d.vel")
      open (8,file="susp_3d.adh")
      open (9,file="susp_3d.orb")

c-----------------------------------------
c read the triangle integration quadrature
c-----------------------------------------

      call gauss_trgl
     +
     +   (mint
     +   ,xiq,etq,wq
     +   )

c----------------------------
c define the collocation nodes
c over the standard triangle
c in the xi-eta plane
c
c compute the inverse
c of the vandermonde matrix
c---------------------------

       call vander (mpoly,Ipd)

c-----------------------
c initialize and prepare
c-----------------------

      Iwall = 0

      if(Iflow.eq.10) Iwall = 1

      Dth = 0.5D0*Dt

      Do i=1,Nprtcl
        theta(i) = theta(i)*pi
          phi(i) =   phi(i)*pi
      End Do

c: global particle interaction velocity

      Do i=1,Nprtcl
       Do k=1,2306
        Do l=1,3
            uclpg(k,l,i) = 0.0D0
             dplg(k,l,i) = 0.0D0
          Do m=1,Nprtcl
            uclpp(k,l,i,m) = 0.0D0
          End Do
         End Do
       End Do
      End Do

c: for adhesion

      Do i=1,2306
        dpl1(i,1) = 0.0D0
        dpl2(i,2) = 0.0D0
        dpl3(i,3) = 0.0D0
        dpl4(i,3) = 0.0D0
      End Do

c: triply periodic flow

      a11 = 1.0D0
      a12 = 0.0D0
      a13 = 0.0D0

      a21 = 0.0D0
      a22 = 1.0D0
      a23 = 0.0D0

      a31 = 0.0D0
      a32 = 0.0D0
      a33 = 1.0D0

c: counters and time

      Kstep = 1           ! global time step counter

      time(1) = 0.0D0

c-.-.-.-.-.-.-.-.-.-.-.-.-.
c time stepping begins here
c-.-.-.-.-.-.-.-.-.-.-.-.-.

 90   Continue

      if(Ienrd.eq.1) then
        read  (2,*)  Nstep
      else if(Ienrd.eq.2) then       ! interactive
        write (6,*)
        write (6,*) 'Enter the number of steps before pausing'
        write (6,*) '0 to quit'
        write (6,*) '---------'
        read  (5,*)  Nstep
      end if

      write (6,*) " susp_3d: Another batch of ",Nstep," steps:"

      if(Nstep.eq.0) Go to 99

c----  another batch:

      Istep  = 1      ! batch step counter

c----  another time step:

  97  Continue

      write (6,*) "------------------------"
      write (6,105) Istep,Nstep,time(Kstep)

c---------------------------------------
c for more than one particle,
c compute the particle collocation nodes
c and store them in a global matrix
c---------------------------------------

       if(Nprtcl.gt.1) then

        Do i=1,Nprtcl   ! loop over particles

         call susp_3d_col
     +
     +   (Ioctaicos
     +   ,ndiv
     +   ,boa(i),coa(i)
     +   ,req(i)
     +   ,Isize
     +   ,cx(i),cy(i),cz(i)
     +   ,theta(i),phi(i)
     +   ,mpoly
     +   ,npts,nelm
     +   ,necl
     +   ,ngcl
     +   ,pcl
     +   )

         Do k=1,ngcl
          Do l=1,3
           pclg(k,l,i) = pcl(k,l)    ! global matrix
          End Do
         End Do

       End Do

       else

        ngcl =  2306

       end if

c--------------------------
c prepare for periodic flow
c--------------------------

      call sgf_3d_3p_ewald
     +
     +  (a11,a12,a13
     +  ,a21,a22,a23
     +  ,a31,a32,a33
     +  ,b11,b12,b13
     +  ,b21,b22,b23
     +  ,b31,b32,b33
     +  ,ew
     +  ,tau
     +  )

      call sgf_3d_3p_vvv
     +
     +   (b11,b12,b13
     +   ,b21,b22,b23
     +   ,b31,b32,b33
     +   ,Max2
     +   ,ew
     +   )

c-----------------------------------
c solve the integral equation 
c and compute the particle translational 
c and angular velocities
c-----------------------------------

      Ipass = 0
 93   Continue        ! outer iterations
      Ipass = Ipass+1

      Do k=1,3   ! initialize the particle stress tensor
       Do l=1,3
         stress(k,l) = 0.0D0
       End Do
      End Do

c-----------------------------
      Do i=1,Nprtcl   ! loop over particles: goulis
c-----------------------------


        Do k=1,ngcl  ! loop over collocation nodes
         Do l=1,3
              uclp(k,l) = uclpg(k,l,i) ! mutually induced particle velocity
               dpl(k,l) =  dplg(k,l,i)  ! initial guess for the dipole
         End Do
c         write (6,104) k,uclp(k,1),uclp(k,2),uclp(k,3)
        End Do

c---
c     write (6,*) " Self-interaction of particle :",i
c---

        call susp_3d_slv
     +
     +   (Ioctaicos
     +   ,ndiv
     +   ,boa(i),coa(i)
     +   ,req(i)
     +   ,Isize
     +   ,cx(i),cy(i),cz(i)
     +   ,theta(i),phi(i)
     +
     +   ,mpoly
     +   ,intm        ! integration method
     +   ,mint
     +   ,tolin
     +   ,Nter
     +
     +   ,uclp      ! particle disturbance velocity
     +
     +   ,Iadhere
     +   ,torque,adhesion
     + 
     +   ,npts,nelm   ! output 
     +   ,p,n
     +   ,area,vlm   !       
     +   ,necl,ngcl
     +   ,pcl            ! collocation points
     +   ,dpl            ! dipole
     +   ,dpldn(i)
     +   ,stresslet
     +   ,Vx(i),Vy(i),Vz(i)
     +   ,Ox(i),Oy(i),Oz(i) 
     +
     +   ,xfadh(i),yfadh(i),zfadh(i)
     +   )

c---
c save the collocation points in a global matrix
c for graphics display
c---

       Do k=1,npts          ! graphics 
        Do l=1,3     
          pg(k,l,i) = p(k,l)
        End Do
       End Do

c---
c update the initial guess (final solution) for the dipole
c---

        Do k=1,ngcl
         Do l=1,3
            dplg(k,l,i)=dpl(k,l)
         End Do
        End Do

c---
c contribution to the particle stress tensor
c---

       Do k=1,3
        Do l=1,3
         stress(k,l) = stress(k,l)+stresslet(k,l)
        End Do
       End Do

c      Do is=1,3
c       write (6,199) (stresslet(is,js),js=1,3)
c      End Do

c--------
c display
c--------

c     area_red = area/(pi4*req(i)**2)          ! normalize
c     vlm_red  = 3.0D0*vlm /(pi4*req(i)**3)    ! normalize
c     write (6,110) area_red
c     write (6,111) vlm_red
c     write (6,800) cx(i),cy(i),cz(i)

c     write (6,*)
c     write (6,*) " Converged solution: x-position, dipole"
c     write (6,*)
c     write (6,*) ngsp
c     Do i=1,ngsp
c      write (6,100) i,pcl(i,1),dpl(i,1),dpl(i,2),dpl(i,3)
c     End Do

      if(Nprtcl.eq.1) Go to 94

c-------------------
c mutual interaction
c-------------------

c----
c: loop over all other particles
c----

       Do j=1,Nprtcl

        If(j.ne.i) then

c       write (6,*) " mutual-interaction of particles :",i,j

        ! set the evaluation points equal to
        ! the jth particle collocation points:

        Do k=1,ngcl 
          Do l=1,3 
             pev(k,l) = pclg(k,l,j)
          End Do
        End Do

        call susp_3d_dlp_mutual
     +
     +   (npts
     +   ,nelm
     +   ,intm        ! integration method
     +   ,mint
     +   ,mpoly
     +   ,necl,ngcl
     +   ,pcl
     +   ,dpl
     +
     +   ,ngcl     ! number of evaluation points
     +   ,pev      ! evaluation points
     +   ,uclp
     +   )

        Do k=1,ngcl      ! update the global vector for the 
          Do l=1,3       ! disturbance velocities
             uclpg(k,l,j) = uclpg(k,l,j)+uclp(k,l)-uclpp(k,l,j,i)
             uclpp(k,l,j,i) = uclp(k,l)
          End Do
c         write (6,*) uclpg(k,1,j),uclpg(k,2,j),uclpg(k,3,j)
        End Do

       End If

      End Do

c     if(Iflow.eq.10) Vx(i) = Vx(i)-shrt*(cy-wall)

c-----------
      End Do   ! goulis
c-----------

      if(Ipass.gt.1)  then

       Do i=1,Nprtcl
       error = (Vx_old(i)-Vx(i))**2+(Vy_old(i)-Vy(i))**2
     +        +(Vz_old(i)-Vz(i))**2
     +        +(Ox_old(i)-Ox(i))**2+(Oy_old(i)-Oy(i))**2
     +        +(Oz_old(i)-Oz(i))**2
       End Do
       error = Dsqrt(error/(3.0D0*Nprtcl))
       write (6,463) Ipass,error
       if(error.lt.tolout) Go to 94

      end if

      Do i=1,Nprtcl
       Vx_old(i) = Vx(i)
       Vy_old(i) = Vy(i)
       Vz_old(i) = Vz(i)
       Ox_old(i) = Ox(i)
       Oy_old(i) = Oy(i)
       Oz_old(i) = Oz(i)
      End Do

      Go to 93

 94   Continue

c------
c print
c------

c infinite space

      Do i=1,Nprtcl

      if(Iwall.eq.0) then
       write (6,456) i,dpldn(i),Vx(i),Vy(i),Vz(i),Ox(i),Oy(i),Oz(i)
      end if

c wall bounded flow:

      if(Iwall.eq.1) then
       Oxr = 2.0D0*Ox(i)
       Oyr = 2.0D0*Oy(i)
       Ozr = 2.0D0*Oz(i)
       dd  = cy(i)-wall
       Vxr = Vx(i)/dd
       Vyr = Vy(i)/dd
       Vzr = Vz(i)/dd
       write (6,456) i,dpldn(i),Vxr,Vyr,Vzr,Oxr,Oyr,Ozr
      end if

      End Do

      write (6,*) " susp_3d: particle stress tensor:"

      Do i=1,3
       Do j=1,3
         stress(i,j) = stress(i,j)/Nprtcl
       End Do
      End Do

c---
      if(Iflow.eq.3.or.Iflow.eq.4.or.Iflow.eq.5
     +      .or.Iflow.eq.20) then

       Do i=1,3
        write (6,199) (stress(i,j),j=1,3)
       End Do

      end if
c---

c-----------------------------------------
c rate of evolution of the director angles
c-----------------------------------------

      Do i=1,Nprtcl

      cth = Dcos(theta(i))
      sth = Dsin(theta(i))
      cph = Dcos(phi(i))
      sph = Dsin(phi(i))

      write (6,*) sth

      dirx(i) = cth      ! particle director
      diry(i) = sth*cph  ! particle director
      dirz(i) = sth*sph  ! particle director

      wx = Oy(i)*dirz(i)-Oz(i)*diry(i)
      wy = Oz(i)*dirx(i)-Ox(i)*dirz(i)
      wz = Ox(i)*diry(i)-Oy(i)*dirx(i)

      wt = -sth*wx + cth*cph*wy + cth*sph*wz
      wp = -sph*wy + cph*wz

      dtheta(i) = wt
c       dphi(i) = 0.0D0
        dphi(i) = wp/sth

      End Do

c---------------
c Jeffery angles
c---------------

      Do i=1,Nprtcl
       chi(i) = Dasin(diry(i)/Dsqrt(dirx(i)**2+diry(i)**2))
       If(dirx(i).lt.0) chi(i) =-chi(i)-pi
c      If(chi(i)/pi.lt.-1.49) Go to 99
       fc = (1.0D0-boa(i)**2)/(1.0D0+boa(i)**2) ! jeffery orbits:

        ! for phi=0:
       dthjef(i) =-0.50D0*(1.0D0 -fc*Dcos(2.0D0*theta(i)) )
      End Do

c---------
c printing
c---------

      write (3,107) Kstep,time(Kstep)
     +            ,(-chi(i)/pi,cx(i),cy(i),cz(i)
     +              ,dirx(i),diry(i),dirz(i)
     +              ,theta(i)/pi,phi(i)/pi
     +              ,i=1,Nprtcl)

      write (6,107) Kstep,time(Kstep)
     +            ,(-chi(i)/pi,cx(i),cy(i),cz(i)
     +              ,dirx(i),diry(i),dirz(i)
     +              ,theta(i)/pi,phi(i)/pi
     +              ,i=1,Nprtcl)

      write (4,107) Kstep,time(Kstep)
     +            ,(-chi(i)/pi
     +            ,Vx(i),Vy(i),Vz(i)
     +            ,Ox(i),Oy(i),Oz(i)
     +            ,-dtheta(i),-dphi(i),-dthjef(i)
     +            ,i=1,Nprtcl)

c     write (4,107) Kstep,time(Kstep)
c    +            ,(-theta(i)/pi
c    +            ,-Oz(i),-dthjef(i)
c    +            ,i=1,Nprtcl)

      write (1,104) Kstep,time(Kstep)
     +              ,stress(1,1),stress(1,2),stress(1,3)
     +                          ,stress(2,2),stress(2,3)
     +                                      ,stress(3,3)

c     write (6,104) Kstep,time(Kstep)
c    +              ,stress(1,1),stress(1,2),stress(1,3)
c    +                          ,stress(2,2),stress(2,3)
c    +                                      ,stress(3,3)

      write (9,102) (dirx(i),diry(i),dirz(i),i=1,Nprtcl)

      write (8,154) Kstep,time(Kstep)
     +   ,(xfadh(i),yfadh(i),zfadh(i),i=1,Nprtcl)

c----------------
c For RK2:
c save coordinates and velocities
c---------------

      if(IRK.eq.2) then

       Do i=1,Nprtcl

        cx_save(i) = cx(i)
        cy_save(i) = cy(i)
        cz_save(i) = cz(i)

        theta_save(i) = theta(i)
          phi_save(i) = phi(i)

        Vx_save(i) = Vx(i)
        Vy_save(i) = Vy(i)
        Vz_save(i) = Vz(i)

        Ox_save(i) = Ox(i)
        Oy_save(i) = Oy(i)
        Oz_save(i) = Oz(i)

        dtheta_save(i) = dtheta(i)
        dphi_save(i) = dphi(i)

       End Do

      end if

c--------------------
c advance with Euler:
c--------------------

      Do i=1,Nprtcl

       cx(i)=cx(i)+Vx(i)*Dt
       cy(i)=cy(i)+Vy(i)*Dt
       cz(i)=cz(i)+Vz(i)*Dt

       theta(i) = theta(i) + dtheta(i)*Dt
       phi(i)   = phi(i)   +   dphi(i)*Dt

      End Do

c: periodic flow
   
      if(Iflow.eq.20) then
       a21 = a21+shrt*Dt
       if(a21.gt. 0.5D0*a11) a21 = a21-a11   ! reset
       if(a21.lt.-a11      ) a21 = a21+a11   ! reset
       Do i=1,Nprtcl
        if(cx(i).gt. 1.5D0*a11) cx(i) = cx(i)-a11   ! reset
        if(cx(i).lt.-1.5D0*a11) cx(i) = cx(i)+a11   ! reset
       End Do
      end if

c---------------------------
      if(IRK.eq.1) Go to 333
c---------------------------

      write (6,*) " susp_3d: second RK2 step"

c----
c for more than one particle,
c compute the particle collocation nodes
c and store them in a global matrix
c---

      if(Nprtcl.gt.1) then

        Do i=1,Nprtcl   ! loop over particles

         call susp_3d_col
     +
     +   (Ioctaicos
     +   ,ndiv
     +   ,boa(i),coa(i)
     +   ,req(i)
     +   ,Isize
     +   ,cx(i),cy(i),cz(i)
     +   ,theta(i),phi(i)
     +   ,mpoly
     +   ,npts,nelm
     +   ,necl,ngcl
     +   ,pcl
     +   )

         Do k=1,ngcl
          Do l=1,3
           pclg(k,l,i) = pcl(k,l)    ! global matrix
          End Do
         End Do

       End Do

      end if

c--------------------------
c prepare for periodic flow
c--------------------------

      call sgf_3d_3p_ewald
     +
     +  (a11,a12,a13
     +  ,a21,a22,a23
     +  ,a31,a32,a33
     +  ,b11,b12,b13
     +  ,b21,b22,b23
     +  ,b31,b32,b33
     +  ,ew
     +  ,tau
     +  )

      call sgf_3d_3p_vvv
     +
     +   (b11,b12,b13
     +   ,b21,b22,b23
     +   ,b31,b32,b33
     +   ,Max2
     +   ,ew
     +   )

c-----------------------------------
c solve the integral equation
c to compute the particle translational
c and angular velocities
c-----------------------------------

      Ipass = 0
 83   Continue        ! outer iterations
      Ipass = Ipass+1

c------------------
      Do i=1,Nprtcl   ! loop over particles
c------------------

        Do k=1,ngcl
         Do l=1,3
           uclp(k,l) = uclpg(k,l,i) ! mutually induced particle velocity
            dpl(k,l) = dplg(k,l,i)  ! initial guess for the dipole
         End Do
        End Do

c     write (6,*) " Self-interaction of particle :",i

        call susp_3d_slv
     +
     +   (Ioctaicos
     +   ,ndiv
     +   ,boa(i),coa(i)
     +   ,req(i)
     +   ,Isize
     +   ,cx(i),cy(i),cz(i)
     +   ,theta(i),phi(i)
     +
     +   ,mpoly
     +   ,intm        ! integration method
     +   ,mint
     +   ,tolin
     +   ,Nter
     +
     +   ,uclp      ! particle disturbance velocity
     +
     +   ,Iadhere
     +   ,torque,adhesion
     +
     +   ,npts,nelm   ! output
     +   ,p,n
     +   ,area, vlm   !
     +   ,necl,ngcl
     +   ,pcl            ! collocation points
     +   ,dpl            ! dipole
     +   ,dpldn(i)
     +   ,stresslet
     +   ,Vx(i),Vy(i),Vz(i)
     +   ,Ox(i),Oy(i),Oz(i)
     +
     +   ,xfadh(i),yfadh(i),zfadh(i)
     +   )

c: update the initial guess

        Do k=1,ngcl
         Do l=1,3
            dplg(k,l,i)=dpl(k,l)
         End Do
        End Do

c---
        if(Nprtcl.eq.1) Go to 84
c---

c-------------------
       Do j=1,Nprtcl      ! mutual interaction
c----

        If(j.ne.i) then

c       write (6,*) " mutual-interaction of particles :",i,j

        Do k=1,ngcl
          Do l=1,3
             pev(k,l)=pclg(k,l,j)   ! evaluation points
          End Do
        End Do

        call susp_3d_dlp_mutual
     +
     +   (npts
     +   ,nelm
     +   ,intm        ! integration method
     +   ,mint
     +   ,mpoly
     +   ,necl,ngcl
     +   ,pcl
     +   ,dpl
     +
     +   ,ngcl     ! number of evaluation points
     +   ,pev      ! evaluation points
     +   ,uclp
     +   )

        Do k=1,ngcl      ! update the global vector for the
          Do l=1,3       ! disturbance velocities
            uclpg(k,l,j)   = uclpg(k,l,j)+uclp(k,l)-uclpp(k,l,j,i)
            uclpp(k,l,j,i) = uclp(k,l)
          End Do
        End Do

        end if

       End Do       ! mutual interaction
c-------------------

c-----------
      End Do   ! loop over particles
c-----------

      if(Ipass.gt.1)  then
       Do i=1,Nprtcl
       error = (Vx_old(i)-Vx(i))**2+(Vy_old(i)-Vy(i))**2
     +        +(Vz_old(i)-Vz(i))**2
     +        +(Ox_old(i)-Ox(i))**2+(Oy_old(i)-Oy(i))**2
     +        +(Oz_old(i)-Oz(i))**2
       End Do
       error = Dsqrt(error/(3.0D0*Nprtcl))
       write (6,463) Ipass,error
       If(error.lt.tolout) Go to 84
      end if

      Do i=1,Nprtcl
       Vx_old(i) = Vx(i)
       Vy_old(i) = Vy(i)
       Vz_old(i) = Vz(i)
       Ox_old(i) = Ox(i)
       Oy_old(i) = Oy(i)
       Oz_old(i) = Oz(i)
      End Do

      Go to 83

c: outer iterations converged

 84   Continue

c-----------------------------------------
c rate of evolution of the director angles
c-----------------------------------------

      Do i=1,Nprtcl

       cth = Dcos(theta(i))
       sth = Dsin(theta(i))
       cph = Dcos(phi(i))
       sph = Dsin(phi(i))

       dirx(i) = cth      ! particle director
       diry(i) = sth*cph  ! particle director
       dirz(i) = sth*sph  ! particle director

       wx = Oy(i)*dirz(i)-Oz(i)*diry(i)
       wy = Oz(i)*dirx(i)-Ox(i)*dirz(i)
       wz = Ox(i)*diry(i)-Oy(i)*dirx(i)
       wt = -sth*wx +cth*cph*wy +cth*sph*wz
       wp = -sph*wy +cph*wz

       dtheta(i) = wt
c      dphi(i)   = 0.0D0
       dphi(i)   = wp/sth

      End Do

c-------------
c advance RK2:
c-------------

      Do i=1,Nprtcl

       cx(i) = cx_save(i)+(Vx_save(i)+Vx(i))*Dth
       cy(i) = cy_save(i)+(Vy_save(i)+Vy(i))*Dth
       cz(i) = cz_save(i)+(Vz_save(i)+Vz(i))*Dth

       theta(i) = theta_save(i) + (dtheta_save(i)+dtheta(i))*Dth
       phi(i)   = phi_save(i)   + (dphi_save(i)  +dphi(i)  )*Dth

      End Do

c------------------------
c End of a time step
c
c Reset counters and time
c------------------------

  333 Continue

      Kstep = Kstep+1
      Istep = Istep+1
      time(Kstep) = time(Kstep-1) + Dt

      write (6,*) ngcl

      if(Istep.le.Nstep) Go to 97

      Go to 90     ! return for another batch of steps

c-----
c Done
c-----

 99   Continue

       write (1,888) ndiv,Isize,intm,mint
     +        ,tolin,tolout,Nter,visc
     +        ,wall,Iadhere,torque
     +        ,adhesion,Uinf,shrt
     +        ,Max1,Max2,mpoly,Ipd,IRK,Dt
       write (3,888) ndiv,Isize,intm,mint
     +        ,tolin,tolout,Nter,visc
     +        ,wall,Iadhere,torque
     +        ,adhesion,Uinf,shrt
     +        ,Max1,Max2,mpoly,Ipd,IRK,Dt
       write (4,888) ndiv,Isize,intm,mint
     +        ,tolin,tolout,Nter,visc
     +        ,wall,Iadhere,torque
     +        ,adhesion,Uinf,shrt
     +        ,Max1,Max2,mpoly,Ipd,IRK,Dt
       write (8,888) ndiv,Isize,intm,mint
     +        ,tolin,tolout,Nter,visc
     +        ,wall,Iadhere,torque
     +        ,adhesion,Uinf,shrt
     +        ,Max1,Max2,mpoly,Ipd,IRK,Dt
       write (9,888) ndiv,Isize,intm,mint
     +        ,tolin,tolout,Nter,visc
     +        ,wall,Iadhere,torque
     +        ,adhesion,Uinf,shrt
     +        ,Max1,Max2,mpoly,Ipd,IRK,Dt

      close (1)
      close (3)
      close (4)
      close (8)
      close (9)

c--------------------
c print all triangles
c--------------------

      open (1,file="trgl.net")

      Index = 2

      if(Index.eq.1) then          ! 6-node triangles
        write (1,*) Nseven
        write (1,*) 7*nelm
        write (1,*) nelm
      else if(Index.eq.2) then     ! 3-node triangles
        write (1,*) Nfour
        write (1,*) 16*nelm
        write (1,*) 4*nelm
      end if

      Do k=1,nelm
       call printel (k,Index,dpl)  ! print in file "caps_3d.net"
      End Do

c print collocation nodes:

      write (1,*) ngcl
c     write (6,*) ngcl

      Do i=1,ngcl
        write (1,102) pcl(i,1),pcl(i,2),pcl(i,3)
c       write (6,102) pcl(i,1),pcl(i,2),pcl(i,3)
      End Do

      close (1)

c------------
c Really Done
c------------

  100 Format (1x,i5,10(1x,f10.5))
  101 Format (f15.10)
  102 Format (10(1x,f12.8))
  104 Format (1x,i4,15(1x,f9.6))
  105 Format (" Step ",i4," out of ",i5,"; time :",F15.10)
  107 Format (1x,i4,1x,f10.7,50(1x,f12.7))
  199 Format (3(1x,f10.7))

  110 Format(" Surface Area :",F15.10)
  111 Format(" Volume       :",F15.10)

  154 Format (1x,i4,1x,f9.6,5(1x,f15.8))

  456 Format (" Particle: ",i2,1x," Proj:",f12.8,/
     +       ," Vx: ",f12.8," Vy: ",f12.8," Vz: ",f12.8,/
     +       ," Ox: ",f12.8," Oy: ",f12.8," Oz: ",f12.8)
  459 Format ("Flow across the boundary:",f15.10)
  463 Format (" Outer iteration: ",I3," error: ", f15.10)

  800 Format(" centroid:",3(1x,f10.5))
  888 Format(/,
     +       " ndiv = ",I3,
     +       " Isize= ",I3,
     +       /,
     +       " Intm = ",I3,
     +       " mint = ",I3,
     +       /,
     +       " tolin = ",f15.10,
     +       " tolout = ",f15.10,
     +       " Nter = ",I3,
     +       /,
     +       " visc = ",f15.10,
     +       " wall = ",f15.10,
     +       /,
     +       " Iadhere = ",I3,
     +       " torque = ",f15.10,
     +       " adhesion = ",f15.10,
     +       /,
     +       " Uinf = ",f15.10,
     +       " shrt = ",f15.10,
     +       /,
     +       " Max1 = ",I3,
     +       " Max2 = ",I3,
     +       /,
     +       " mpoly = ",I3,
     +       " Ipd = ",I3,
     +       /,
     +       " IRK = ",I3,
     +       " Dt = ",f15.10
     +       )

      stop
      end
