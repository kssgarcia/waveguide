#
# Objects
# -------
#
OBJ1  = trgl6_octa_ss.o gauss_trgl.o gauss_leg.o
OBJ2  = sgf_3d_w.o
OBJ3  = bump_3d.o bump_3d_slp.o bump_3d_slp_sing.o bump_3d_slp_sing_lin.o
OBJ31 = elm_geom.o interp_p.o abc.o printel.o
OBJ4  = gel_inv.o gel.o
OBJs  = $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ31) $(OBJ4) $(OBJ5) 
#
# link
# ----
#
bump_3d: $(OBJs)
	f77 -o bump_3d $(OBJs)
#
# compile
# ------
#
bump_3d.o: bump_3d.f 
	f77 -c bump_3d.f 
bump_3d_slp.o: bump_3d_slp.f 
	f77 -c bump_3d_slp.f 
bump_3d_slp_sing.o: bump_3d_slp_sing.f 
	f77 -c bump_3d_slp_sing.f 
bump_3d_slp_sing_lin.o: bump_3d_slp_sing_lin.f 
	f77 -c bump_3d_slp_sing_lin.f 
gauss_leg.o: gauss_leg.f 
	f77 -c gauss_leg.f 
gauss_trgl.o: gauss_trgl.f 
	f77 -c gauss_trgl.f 
trgl6_octa_ss.o: trgl6_octa_ss.f 
	f77 -c trgl6_octa_ss.f 
sgf_3d_w.o: sgf_3d_w.f 
	f77 -c sgf_3d_w.f 
gel_inv.o: gel_inv.f 
	f77 -c gel_inv.f 
gel.o: gel.f 
	f77 -c gel.f 
elm_geom.o: elm_geom.f 
	f77 -c elm_geom.f 
interp_p.o: interp_p.f
	f77 -c interp_p.f
abc.o: abc.f
	f77 -c abc.f
printel.o: printel.f
	f77 -c printel.f
#
# all
# ---
#
all:
	make bump_3d
#
# clean
# -----
#
clean:
	rm -f core PLOTDAT
	rm -f $(OBJs) bump_3d
	rm -f particle.points
	rm -f bump_3d.net
#
# purge
# -----
#
purge:
	rm -f core PLOTDAT
	rm -f bump_3d.net
	rm -f $(OBJs)
	rm -f particle.points
