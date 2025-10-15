      subroutine prtcl_ax_str
     +
     +  (Iflow
     +  ,X00,Y00
     +  ,Mstr
     +  ,IRK
     +  ,Dl
     +  ,L
     +  )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------
c Compute one streamline
c starting at the point: X00,Y00
c
c Legend:
c ------
c
c Mstr: Maximum number of points
c
c Capacity:
c --------
c
c     25 particles
c     64 elements along each particle
c------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension NE(25),Itp(25)
      Dimension xw(25,65),yw(25,65),tw(25,65)

      Dimension aux1(65),aux2(65),aux3(65)

      Dimension xstr(900),ystr(900) ! for streamlines

c--------------
c common blocks
c--------------

      common/INTGR1/Nprtcl,NGL,NE,Itp,Ncl

      common/points/xw,yw,tw
      common/stream/xstr,ystr

c---------------------------------
c Check whether the starting point
c lies inside a particle
c by the circulation method
c---------------------------------

      Do i=1,Nprtcl

       Do j=1,NE(i)+1
        aux1(j) = xw(i,j)-X00
        aux2(j) = yw(i,j)-Y00
        aux3(j) = Dsqrt(aux1(j)**2+aux2(j)**2)
       End Do

       angle = 0.0D0

       Do j=1,NE(i)

        If(aux3(j).gt.0) then
          tmp = aux1(j)*aux1(j+1)+aux2(j)*aux2(j+1)
          tmp = tmp/(aux3(j)*aux3(j+1))
          dangle = acos(tmp)
          tmp = aux1(j)*aux2(j+1)-aux2(j)*aux1(j+1)
          If(tmp.lt.0) dangle = - dangle
          angle = angle+dangle
        End If

       End Do

c-----------------------------------------
c treat singly-connected particle contours
c with the points 1 and NE+1
c lying on the x axis
c----------------------------------------

       If(abs(aux3(1)-aux3(NE(i)+1)).gt.0.000001) then

        j = NE(i)+2
        xw(i,j) = xw(i,1)
        yw(i,j) = yw(i,1)
        aux1(j) = xw(i,j)-X00
        aux2(j) = yw(i,j)-Y00
        aux3(j) = sqrt(aux1(j)**2+aux2(j)**2)
        j = NE(i) + 1
        tmp = aux1(j)*aux1(j+1)+aux2(j)*aux2(j+1)
        tmp = tmp/(aux3(j)*aux3(j+1))
        dangle = acos(tmp)
        tmp = aux1(j)*aux2(j+1)-aux2(j)*aux1(j+1)
        If(tmp.lt.0) dangle = - dangle
        angle = angle+dangle
       End If

c---
c check the circulation
c---

       If(abs(angle).gt.0.10) then
         write (6,*) " prtcl_ax_str: starting point inside a particle"
         L = 0
         Go to 99
       End If

      End Do

c-------------------------------
c End of checking
c
c Proceed to compute streamlines
c-------------------------------

      xstr(1) = X00
      ystr(1) = Y00

      K = 1     ! local counter for inquiry
      L = 1     ! total counter

  20  Continue

      xstr(L) = X00
      ystr(L) = Y00

c---
c first velocity evaluation
c---

      call prtcl_ax_vel
     +
     +   (Iflow
     +   ,X00,Y00
     +   ,Ux1,Uy1
     +   )

      write (6,104) L,X00,Y00,Ux1,Uy1

      step = Dl/dsqrt(Ux1**2+Uy1**2) ! time step

      Xsv = X00  ! save
      Ysv = Y00

      steph = 0.5D0*step

c-----------------
      If(IRK.eq.2) then
c-----------------

        X00 = Xsv + step * Ux1
        Y00 = Ysv + step * Uy1

        call prtcl_ax_vel
     +
     +    (Iflow
     +    ,X00,Y00
     +    ,Ux2,Uy2
     +    )

        X00 = Xsv + steph*(Ux1+Ux2)
        Y00 = Ysv + steph*(Uy1+Uy2)
c----------------------
      Else If(IRK.eq.4) then
c----------------------

        X00 = Xsv + steph * Ux1
        Y00 = Ysv + steph * Uy1

        call prtcl_ax_vel
     +
     +     (Iflow
     +     ,X00,Y00
     +     ,Ux2,Uy2
     +     )

        X00 = Xsv + steph * Ux2
        Y00 = Ysv + steph * Uy2

        call prtcl_ax_vel
     +
     +      (Iflow
     +      ,X00,Y00
     +      ,Ux3,Uy3
     +      )

        X00 = Xsv + step * Ux3
        Y00 = Ysv + step * Uy3

        call prtcl_ax_vel
     +
     +      (Iflow
     +      ,X00,Y00
     +      ,Ux4,Uy4
     +      )

        X00 = Xsv + step * (Ux1 +2.0D0*Ux2 +2.0D0*Ux3 +Ux4)/6.0D0
        Y00 = Ysv + step * (Uy1 +2.0D0*Uy2 +2.0D0*Uy3 +Uy4)/6.0D0
c---
      End If
c---

      K = K+1
      L = L+1

      xstr(L) = X00
      ystr(L) = Y00

c----------------
c Stopping checks
c----------------

c----------------------------------
c Check for a closed streamline
c by examining the distance between the
c current and the first point
c----------------------------------

      If(L.gt.2) then
       test = Dsqrt((X00-Xstr(1))**2+(Y00-Ystr(1))**2)
       if(test.lt.stpsp) Go to 99
      End If

c---------------------------
      If(K.lt.Mstr) Go to 20
c---------------------------

      K = 1    ! reset local counter

      write (6,*) 
      write (6,*) " Continue this streamline ?"
      write (6,*) 
      write (6,*) " Enter 0 for No, 1 for yes"
      write (6,*) " -------------------------"

      read  (5,*) Icon

      If(Icon.eq.1) Go to 20

c-----
c Done
c-----

  99  Continue

 104  Format (1x,i3,20(1x,f9.5))

      Return
      End
