      program fin

      Implicit Double Precision (a-h,o-z)

      Dimension z(0:500,0:500),r(0:500,0:500)
      Dimension T(0:500,0:500)
      Dimension Tnew(0:500,0:500)
      Dimension TI1(0:500),TI2(0:500)

      Double Precision k1,k2

      Parameter (Nsteps=1000000)

c-------------------------------------------
c Evolution of temperature distribution
c in a fin.
c
c SYMBOLS:
c -------
c
c k1: thermal conductivity of aluminum
c k2: thermal conductivity of silica
c
c hc: thermal conductance of the interface
c-------------------------------------------

c----------
c constants
c----------

      Null = 0
      Ntwo = 2

c-----------------------------
c all values in the MKS system
c-----------------------------

      a      = 0.05
      b      = 0.001
      k1     = 177.0
      rho1   = 2707.0
      cp1    = 892.00
      k2     = 1.38
      rho2   = 1332.00
      cp2    = 745.00
      hc     = 100.00
      Tw     = 250.00
      Tinit  = 50.00

c----------
c grid size
c----------

      M = 4*4
      N = 4*4

c------------------------------
c set the numerical diffusivity
c------------------------------

      alpha1 = 0.3

c---------------
c inquire input
c---------------
c
c     write (6,*)
c     write (6,*) " Please enter a"
c     write (6,*) " --------------"
c     read  (5,*) a

c     write (6,*)
c     write (6,*) " Please enter b"
c     write (6,*) " --------------"
c     read  (5,*) b

c     write (6,*)
c     write (6,*) " Please enter k1"
c     write (6,*) " ---------------"
c     read  (5,*) k1

c     write (6,*)
c     write (6,*) " Please enter rho1"
c     write (6,*) " ---------------"
c     read  (5,*) rho1

c     write (6,*)
c     write (6,*) " Please enter cp1"
c     write (6,*) " ---------------"
c     read  (5,*) cp1

c     write (6,*)
c     write (6,*) " Please enter k2"
c     write (6,*) " ---------------"
c     read  (5,*) k2

c     write (6,*)
c     write (6,*) " Please enter rho2"
c     write (6,*) " ---------------"
c     read  (5,*) rho2

c     write (6,*)
c     write (6,*) " Please enter cp2"
c     write (6,*) " ---------------"
c     read  (5,*) cp2

c     write (6,*)
c     write (6,*) " Please enter the interface thermal conductance hc"
c     write (6,*) " -------------------------------------------------"
c     read  (5,*) hc

c     write (6,*)
c     write (6,*) " Please enter Tw"
c     write (6,*) " ---------------"
c     read  (5,*) Tw

c     write (6,*)
c     write (6,*) " Please enter M"
c     write (6,*) " --------------"
c     read  (5,*) M

c     write (6,*)
c     write (6,*) " Please enter N"
c     write (6,*) " --------------"
c     read  (5,*) N

c     write (6,*)
c     write (6,*) " Please enter alpha1"
c     write (6,*) " -------------------"
c     read  (5,*) alpha1

c--------
c prepare
c--------

      Dz = a/M
      Dr = b/N

      M1 = M+1
      M2 = M+2

      N1 = N+1
      N2 = N+2
      N4 = N/4

      rat = k2/k1

      beta = (Dr/Dz)**2

c----------
c time step
c----------

      Dt = alpha1*rho1*cp1*Dr**2/k1

      alpha2 = Dt*k2/(rho2*cp2*Dr**2)

c-------------------------
c generate the grid points
c-------------------------

      Do i=0,M2
       zi = (i-1.0)*Dz
       Do j=0,N2
        z(i,j) = zi
        r(i,j) = (j-1.0)*Dr
       End Do
      End Do

c-----------
c initialize
c-----------
   
      time = 0.0D0

      Do i=1,M+2
       Do j=1,N+2
        T(i,j) = Tinit
       End Do
      End Do

      Do i=1,M+2
        TI1(i) = Tinit
        TI2(i) = Tinit
      End Do

      open (4,file="fin.out1")

c-----------
c time stepping
c-----------

       Icount = 1

       write (4,101) Icount,time,T(M1,1)
       write (6,101) Icount,time,T(M1,1)

       Jcount = 1
       write (4,101) Jcount,time,T(M1,1)
       write (6,101) Jcount,time,T(M1,1)

       Ipoint = 0

  97   Continue

C-----------------------
       Do Istep=1,Nsteps
C-----------------------

c------
c left wall
c------

       Do j=0,N+2
        T(1,j) = Tw
       End Do

c------
c upper wall
c------

       Do i=1,M+2
        T(i,N+2) = T(i,N)
       End Do

c------
c lower wall
c------

       Do i=1,M+2
        T(i,0) = T(i,2)
       End Do

c------
c right wall
c------

       Do j=0,N+2
        T(M+2,j) = T(M,j)
       End Do

c----
c interface
c-----

       Do i=1,M+2

       a11 = 1.0D0
       a12 = rat
       b1  = T(i,N4)+ rat * T(i,N4+2)

       a21 =  k1/Dr + hc
       a22 = -hc
       b2  = k1/Dr * T(i,N4)

       Det = a11*a22-a12*a21

       TI1(i) = (b1*a22-b2*a12)/Det
       TI2(i) = (b2*a11-b1*a21)/Det

       End Do

c----
c Differential equation
c-----
    
        Do i=2,M1
         j = 1
         Tnew(i,j) = T(i,j)
     +          + alpha1*(
     +               +2.0*(T(i,j+1)-2.0*T(i,j)+T(i,j-1))
     +         +beta*(T(i+1,j)-2.0*T(i,j)+T(i-1,j))
     +               ) 
         Do j=2,N4-1
          Tnew(i,j) = T(i,j)
     +           + alpha1*(
     +                0.5*Dr/r(i,j) *(T(i,j+1)-T(i,j-1))
     +                +T(i,j+1)-2.0*T(i,j)+T(i,j-1)
     +          +beta*(T(i+1,j)-2.0*T(i,j)+T(i-1,j))
     +                ) 
         End Do
         j = N4
         Tnew(i,j) = T(i,j)
     +           + alpha1*(
     +                0.5*Dr/r(i,j) *(TI1(i)-T(i,j-1))
     +                +TI1(i)-2.0*T(i,j)+T(i,j-1)
     +          +beta*(T(i+1,j)-2.0*T(i,j)+T(i-1,j))
     +                ) 
         j = N4+2
         Tnew(i,j) = T(i,j)
     +           + alpha2*(
     +                0.5*Dr/r(i,j) *(T(i,j+1)-TI2(i))
     +                +T(i,j+1)-2.0*T(i,j)+TI2(i)
     +          +beta*(T(i+1,j)-2.0*T(i,j)+T(i-1,j))
     +                ) 
         Do j=N4+3,N+1
          Tnew(i,j) = T(i,j)
     +           + alpha2*(
     +                0.5*Dr/r(i,j) *(T(i,j+1)-T(i,j-1))
     +                +T(i,j+1)-2.0*T(i,j)+T(i,j-1)
     +          +beta*(T(i+1,j)-2.0*T(i,j)+T(i-1,j))
     +                ) 
         End Do

        End Do

c-------
c update
c-------

       Do i=2,M1
        Do j=1,N1
         T(i,j) = Tnew(i,j)
        End Do
       End Do

       time = time+Dt

       Icount = Icount + 1
       Jcount = Jcount + 1

c      Tmin = 1000.0
c      Do i=1,M1
c       If(TI1(i).lt.Tmin) Tmin = TI1(i)
c      End Do

       If(Jcount.eq.5000) then
         write (4,101) Ipoint,time,T(M1,N1)
         write (6,101) Ipoint,time,T(M1,N1)
c        write (4,101) Ipoint,time,Tmin
c        write (6,101) Ipoint,time,Tmin
         Jcount = 1
         Ipoint = Ipoint + 1
       End If

C------------
       End Do
C------------

       write (6,*)
       write (6,*) " Continue the time stepping?"
       write (6,*) " Enter 0 for NO"
       read  (5,*) Icon

       If(Icon.ne.0) Go to 97
       
c----------
c print out
c----------

       open (9,file="fin.out2")

c---------------
c print the grid
c---------------

       Go to 333

       Do i=1,M1
        write (9,*) Ntwo
        write (9,100) z(i,1) /a,r(i,1)/b,Null
        write (9,100) z(i,N1)/a,r(i,N1)/b,Null
       End Do

       Do j=1,N1
        write (9,*) Ntwo
        write (9,100) z(1,j)/a,r(1,j)/b,Null
        write (9,100) z(M1,j)/a,r(M1,j)/b,Null
       End Do

 333   Continue

c------------
c temperature
c------------

       Do i=1,M1
        write (9,*) N2
        Do j=1,N4
         write (9,100) z(i,j)/a,r(i,j)/b,T(i,j)/Tw
        End Do
         j = N4+1
         write (9,100) z(i,j)/a,r(i,j)/b,TI1(i)/Tw
         write (9,100) z(i,j)/a,r(i,j)/b,TI2(i)/Tw
        Do j=N4+2,N+1
         write (9,100) z(i,j)/a,r(i,j)/b,T(i,j)/Tw
        End Do
       End Do

       Do j=1,N4
        write (9,*) M1
        Do i=1,M+1
         write (9,100) z(i,j)/a,r(i,j)/b,T(i,j)/Tw
        End Do
       End Do

       j = N4+1
       write (9,*) M1
       Do i=1,M1
         write (9,100) z(i,j)/a,r(i,j)/b,TI1(i)/Tw
       End Do
       j = N4+1
       write (9,*) M1
       Do i=1,M1
         write (9,100) z(i,j)/a,r(i,j)/b,TI2(i)/Tw
       End Do

       Do j=N4+2,N1
        write (9,*) M1
        Do i=1,M1
         write (9,100) z(i,j)/a,r(i,j)/b,T(i,j)/Tw
        End Do
       End Do

       write (9,*) Null
       close (9)

c-----
c Done
c-----

  100  Format (3(1x,f10.5))
  101  Format (1x,i9,3(1x,f10.5))

       Stop
       End
