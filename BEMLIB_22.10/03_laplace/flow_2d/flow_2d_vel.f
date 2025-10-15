      subroutine velocity 
     +
     +   (X00,Y00
     +   ,Ux,Uy
     +   )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c--------------------------------------------
c Computes velocity from the bounday-integral
c representation
c
c Velocity is computed by numerical differentiation
c of the potential
c--------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension NE(10),Itp(10)
      Dimension actis(10),xcntr(10),ycntr(10)

      Dimension xw(10,200),yw(10,200),tw(10,200)
      Dimension phi(10,200)

      Dimension dphidn0(500)

      Parameter (eps=0.01D0)

c--------------
c common blocks
c--------------

      common/xxx01/Iflow,NSG,NGL,NE,Itp
      common/xxx02/xw,yw,tw,phi,dphidn0
      common/xxx03/actis,xcntr,ycntr
      common/xxx04/Vx,Vy

c----------
c constants
c----------

      eps2 = 2.0D0*eps

      Ising = 0

c-----------
c initialize
c-----------

      phi1 = 0.0D0
      phi2 = 0.0D0
      phi3 = 0.0D0
      phi4 = 0.0D0

      unused = 0.0D0

      X01 = X00-eps
      X02 = X00+eps

      Y01 = Y00-eps
      Y02 = Y00+eps

c---------------------------------
c Boundary integral representation
c---------------------------------

      j = 0         ! counter

      Do K=1,NSG

       rad  = actis(k)
       xcnt = xcntr(k)
       ycnt = ycntr(k)

        Do L=1,NE(K)

        X1 = XW(K,L)
        Y1 = YW(K,L)
        T1 = TW(K,L)

        X2 = XW(K,L+1)
        Y2 = YW(K,L+1)
        T2 = TW(K,L+1)

        j = j+1

         call flow_2d_sdlp
     +
     +       (Iflow
     +       ,X01,Y00,unused
     +       ,X1,Y1,T1
     +       ,X2,Y2,T2
     +       ,NGL
     +       ,Ising
     +       ,Itp(k)
     +       ,rad,xcnt,ycnt
     +       ,QQQ
     +       ,WWW
     +       )

        phi1 = phi1 - dphidn0(j)*QQQ + phi(k,l)*WWW

         call flow_2d_sdlp
     +
     +       (Iflow
     +       ,X02,Y00,unused
     +       ,X1,Y1,T1
     +       ,X2,Y2,T2
     +       ,NGL
     +       ,Ising
     +       ,Itp(k)
     +       ,rad,xcnt,ycnt
     +       ,QQQ
     +       ,WWW
     +       )

        phi2 = phi2 - dphidn0(j)*QQQ + phi(k,l)*WWW

         call flow_2d_sdlp
     +
     +       (Iflow
     +       ,X00,Y01,unused
     +       ,X1,Y1,T1
     +       ,X2,Y2,T2
     +       ,NGL
     +       ,Ising
     +       ,Itp(k)
     +       ,rad,xcnt,ycnt
     +       ,QQQ
     +       ,WWW
     +       )

        phi3 = phi3 - dphidn0(j)*QQQ + phi(k,l)*WWW

         call flow_2d_sdlp
     +
     +      (Iflow
     +      ,X00,Y02,unused
     +      ,X1,Y1,T1
     +      ,X2,Y2,T2
     +      ,NGL
     +      ,Ising
     +      ,Itp(k)
     +      ,rad,xcnt,ycnt
     +      ,QQQ
     +      ,WWW
     +      )

        phi4 = phi4 - dphidn0(j)*QQQ + phi(k,l)*WWW

        End Do

      End Do

c-------
c numerical differentiation
c-------

      Ux = (phi2-phi1)/eps2 
      Uy = (phi4-phi3)/eps2

c-------------------------
c Add the unperturbed flow
c-------------------------

      Ux = Ux+Vx
      Uy = Uy+Vy

c-----
c Done
c-----

      Return
      End
