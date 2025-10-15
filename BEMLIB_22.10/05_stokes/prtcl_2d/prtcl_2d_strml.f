      subroutine prtcl_2d_strml
     +
     +   (Iflow
     +   ,X0,Y0
     +   ,Mstr,IRK
     +   ,Dl
     +   ,L
     +   ,xstr,ystr
     +   )

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-----------------------------------
c Computation of a streamline
c originating from the point: X0,Y0
c
c Capacity:
c --------
c
c 72 particles
c 64 collocation points
c 900 points along each streamline
c-----------------------------------
 
      Implicit Double Precision (a-h,o-z)

      Dimension NE(72),Itp(72)

      Dimension xw(72,65),yw(72,65),tw(72,65)

      Dimension aux1(65),aux2(65),aux3(65)

      Dimension xstr(900),ystr(900)        ! for streamlines

c--------------
c common blocks
c--------------

      common/INTGR1/Nprtcl,NE,Itp,Ncl
      common/points/xw,yw,tw

      common/aaaa/a11,a12,a21,a22

c---------------------------------
c check whether the starting point 
c lies inside a particle
c by the method of circulation
c---------------------------------

      Do i=1,Nprtcl

       Do j=1,NE(i)+1
        aux1(j) = xw(i,j)-X0
        aux2(j) = yw(i,j)-Y0
        aux3(j) = sqrt(aux1(j)**2+aux2(j)**2)
       End Do

       angle = 0.0D0

       Do j=1,NE(i)
        tmp = aux1(j)*aux1(j+1)+aux2(j)*aux2(j+1)
        tmp = tmp/(aux3(j)*aux3(j+1))
        Dangle = acos(tmp)
        tmp = aux1(j)*aux2(j+1)-aux2(j)*aux1(j+1)
        If(tmp.lt.0) Dangle = - Dangle
        angle = angle+dangle

       End Do

       If(abs(angle).gt.0.10) then
         write (6,*)
         write (6,*) " streamline: One starting point inside a particle"
         L = 0
         Go to 99
       End If

c--------------------------
c For doubly-periodic flow,
c the particle replicas
c--------------------------

       If(Iflow.eq.10) then

       Do j=1,NE(i)+1
        aux1(j) = xw(i,j)-X0+a11
        aux2(j) = yw(i,j)-Y0
        aux3(j) = sqrt(aux1(j)**2+aux2(j)**2)
       End Do

       angle = 0.0

       Do j=1,NE(i)
        tmp = aux1(j)*aux1(j+1)+aux2(j)*aux2(j+1)
        tmp = tmp/(aux3(j)*aux3(j+1))
        Dangle = acos(tmp)
        tmp = aux1(j)*aux2(j+1)-aux2(j)*aux1(j+1)
        If(tmp.lt.0) Dangle = - Dangle
        angle = angle+dangle
       End Do

       If(abs(angle).gt.0.10) then
         write (6,*)
         write (6,*) " One starting point inside a particle"
         L = 0
         Go to 99
       End If

       Do j=1,NE(i)+1
        aux1(j) = xw(i,j)-X0
        aux2(j) = yw(i,j)-Y0+a22
        aux3(j) = sqrt(aux1(j)**2+aux2(j)**2)
       End Do

       angle = 0.0

       Do j=1,NE(i)
        tmp = aux1(j)*aux1(j+1)+aux2(j)*aux2(j+1)
        tmp = tmp/(aux3(j)*aux3(j+1))
        Dangle = acos(tmp)
        tmp = aux1(j)*aux2(j+1)-aux2(j)*aux1(j+1)
        If(tmp.lt.0) Dangle = - Dangle
        angle = angle+dangle
       End Do

       If(abs(angle).gt.0.10) then
         write (6,*) " streamlines: One starting point"
         write (6,*) "              inside a particle"
         L = 0
         Go to 99
       End If

       End If

      End Do

c-------------------------------
c End of checking
c
c Proceed to compute streamlines
c-------------------------------

      xstr(1) = X0
      ystr(1) = Y0

      X = X0
      Y = Y0

      K = 1     ! local counter for inquiry
      L = 1     ! total counter

  20  Continue

c---
c first velocity evaluation
c---

      Xsv = X  ! save
      Ysv = Y

      call velocity 
     +
     +   (Iflow
     +   ,X,Y
     +   ,Ux1,Uy1
     +   )

      step  = Dl/dsqrt(Ux1*Ux1+Uy1*Uy1) ! time step
      steph = 0.5D0*step

c----------------------
      If(IRK.eq.2) then
c----------------------

        X = Xsv + step * Ux1
        Y = Ysv + step * Uy1

        call velocity 
     +
     +    (Iflow
     +    ,X,Y
     +    ,Ux2,Uy2
     +    )

        X = Xsv + steph*(Ux1+Ux2)
        Y = Ysv + steph*(Uy1+Uy2)

c---------------------------
      Else If(IRK.eq.4) then
c---------------------------

        X = Xsv + steph * Ux1
        Y = Ysv + steph * Uy1

        call velocity
     +
     +     (Iflow
     +     ,X,Y
     +     ,Ux2,Uy2
     +     )

        X = Xsv + steph * Ux2
        Y = Ysv + steph * Uy2

        call velocity
     +
     +     (Iflow
     +     ,X,Y
     +     ,Ux3,Uy3
     +     )

        X = Xsv + step * Ux3
        Y = Ysv + step * Uy3

        call velocity 
     +
     +     (Iflow
     +     ,X,Y
     +     ,Ux4,Uy4
     +     )

        X = Xsv + step*(Ux1 +2.0D0*Ux2 +2.0D0*Ux3 +Ux4)/6.0D0
        Y = Ysv + step*(Uy1 +2.0D0*Uy2 +2.0D0*Uy3 +Uy4)/6.0D0

c-----------
      End If
c-----------

      K = K+1
      L = L+1

c     write (6,104) L,X,Y

      xstr(L) = X
      ystr(L) = Y

c-------------------------
c  Perform stopping checks
c-------------------------

c----------------------------------
c Check whether the streamline has closed
c upon itself,
c by examining the distance between the
c current point and the first point
c----------------------------------

      If(L.gt.2) then

       test = dsqrt((X-X0)**2+(Y-Y0)**2)

       If(test.lt.Dl) Go to 99

      End If

c----------------------------------
c In the case of doubly-periodic flow,
c check whether the streamline has
c returned to the original x or y 
c position
c----------------------------------

      If(Iflow.eq.10) then

       testx = abs(X-X0)

       If(testx.gt.a11) then
         xstr(L) = xstr(1)+a11
         ystr(L) = ystr(1)
         Go to 99
       End If

       testy = dabs(Y0-Y0)

       If(testy.gt.a22) then
         xstr(L) = xstr(1)
         ystr(L) = ystr(1)+a22
         Go to 99
       End If

      End If

c-----------------------
c End of stopping checks
c-----------------------

      If(K.lt.Mstr) Go to 20

      K = 1    ! reset local counter

      write (6,*) 
      write (6,*) " Continue this streamline ?"
      write (6,*) 
      write (6,*) " Enter 0 for No, 1 for yes"
      write (6,*) "--------------------------"

      read  (5,*) Icon

      If(Icon.eq.1) Go to 20

c-----
c Done
c-----

  99  Continue

 104  Format (1x,i3,20(1x,f9.5))

      Return
      End
