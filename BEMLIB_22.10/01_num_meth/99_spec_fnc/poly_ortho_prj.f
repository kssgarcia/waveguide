      function poly_ortho_prj (menu,B,i,j)

c=========================================
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c------------------------------------------------
c calculate the projection of x^i onto
c the jth orthogonal polynomial
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension B(0:21,0:21),p(0:21)

      Do k=0,j
        if(menu.eq.1) p(k) = poly_ortho_leg(i+k)
        if(menu.eq.2) p(k) = poly_ortho_rad(i+k)
        if(menu.eq.3) p(k) = poly_ortho_log(i+k)
        if(menu.eq.4) p(k) = poly_ortho_lin(i+k)
        if(menu.eq.5) p(k) = poly_ortho_lgr(i+k)
        if(menu.eq.6) p(k) = poly_ortho_her(i+k)
        Do l=0,k-1
          p(k) = p(k) + B(k,l)*p(l)
        End Do
      End Do

      poly_ortho_prj = p(j)

c-----
c done
c-----

      return
      end
