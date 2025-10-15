      subroutine sm_ax
     +
     +  (NSG
     +  ,P
     +  )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-----------------------------
c Smoothing of property P
c using the 5-point formula of
c Longuett-Higgins and Cokelet
c-----------------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision P(0:900),Ps(0:902)

c--------
c prepare
c--------

      NSGa = NSG-1
      NSG1 = NSG+1
      NSG2 = NSG+2
      NSG3 = NSG+3

c----------------
c Dump and extend
c----------------

      Do i=1,NSG1
       Ps(i) = P(i)
      End Do

      Ps(0)    = Ps(2)
      Ps(NSG2) = Ps(NSG)
      Ps(NSG3) = Ps(NSGa)

c----------
c smoothing
c----------

      P(1) = - Ps(3)+4.0*Ps(2)+10.0*Ps(1) 
     +              +4.0*Ps(2)     -Ps(3)

      Do i=2,NSG1
       P(i) = - Ps(i-2)+4.0*Ps(i-1)+10.0*Ps(i) 
     +                 +4.0*Ps(i+1)     -Ps(i+2)
      End Do

      Do i=1,NSG1
       P(i) = P(i)/16.00
      End Do

c---------
c printing
c---------

      Do i=1,NSG1
       write (6,100) i,Ps(i),P(i)
      End Do

c-----
c Done
c-----

  100 Format (1x,i5,2(1x,f15.10))

      Return
      End
