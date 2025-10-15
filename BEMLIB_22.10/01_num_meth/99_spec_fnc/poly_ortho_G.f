      function poly_ortho_G (menu,B,G,n)

c=========================================
c Copyright by C. Pozrikidis, 1999
c
c All rights reserved
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c----------------------------------------------------
c calculate the norm of the orthogonal polynomial pn
c Gn = <pn|pn>
c----------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension p(0:21)
      Dimension B(0:21,0:21),G(0:21)

c----------------------
c project <x**n | x**n>
c----------------------

      if(menu.eq.1) Gn = poly_ortho_leg(2*n)
      if(menu.eq.2) Gn = poly_ortho_rad(2*n)
      if(menu.eq.3) Gn = poly_ortho_log(2*n)
      if(menu.eq.4) Gn = poly_ortho_lin(2*n)
      if(menu.eq.5) Gn = poly_ortho_lgr(2*n)
      if(menu.eq.6) Gn = poly_ortho_her(2*n)

c---------------------------------------
c Project x**n onto p(k), where k<n
c calculate the cross terms <x**n | pk>
c---------------------------------------

      Do k=0,n-1
          if(menu.eq.1) p(k) = poly_ortho_leg(n+k)
          if(menu.eq.2) p(k) = poly_ortho_rad(n+k)
          if(menu.eq.3) p(k) = poly_ortho_log(n+k)
          if(menu.eq.4) p(k) = poly_ortho_lin(n+k)
          if(menu.eq.5) p(k) = poly_ortho_lgr(n+k)
          if(menu.eq.6) p(k) = poly_ortho_her(n+k)
        Do l=0,k-1
          p(k) = p(k) + B(k,l)*p(l)
        End Do
        Gn = Gn + 2.0D0*B(n,k)*p(k)
      End Do

c--------------------
c calculate <pn | pn>
c--------------------

      Do k=0,n-1
         Gn = Gn + B(n,k)**2 * G(k)
      End Do

      poly_ortho_G = Gn

c-----
c done
c-----

      return
      end
