#
# Objects
# -------
#
OBJ0 = verbal.o
OBJ1 = trgl6_octa.o trgl6_icos.o gauss_trgl.o gauss_leg.o 
OBJ2 = sgf_3d_fs.o sgf_3d_w.o sgf_3d_2p_w.o
OBJ2A = sgf_3d_3p.o sgf_3d_3p_ewald.o sgf_3d_3p_qqq.o
OBJ3 = prtcl_3d_mob.o
OBJ30 = elm_geom.o abc.o interp_p.o printel.o
OBJ33 = slp_trgl6.o slp_trgl6_sing.o slp_trgl3_sing.o
OBJ4 = gel.o gel_inv.o
OBJ  = $(OBJ0) $(OBJ1) $(OBJ2) $(OBJ2A) $(OBJ3) $(OBJ30) $(OBJ33) $(OBJ4)
#
# link
# ----
#
prtcl_3d_mob: $(OBJ)
	f77 -o prtcl_3d_mob $(OBJ)
#
# compile
# ------
#
prtcl_3d_mob.o: prtcl_3d_mob.f 
	f77 -c prtcl_3d_mob.f 
trgl6_octa.o: trgl6_octa.f 
	f77 -c trgl6_octa.f 
trgl6_icos.o: trgl6_icos.f 
	f77 -c trgl6_icos.f 
verbal.o: verbal.f
	f77 -c verbal.f 
sgf_3d_fs.o: sgf_3d_fs.f 
	f77 -c sgf_3d_fs.f 
sgf_3d_w.o: sgf_3d_w.f 
	f77 -c sgf_3d_w.f 
sgf_3d_3p.o: sgf_3d_3p.f 
	f77 -c sgf_3d_3p.f 
sgf_3d_3p_ewald.o: sgf_3d_3p_ewald.f 
	f77 -c sgf_3d_3p_ewald.f 
sgf_3d_3p_qqq.o: sgf_3d_3p_qqq.f 
	f77 -c sgf_3d_3p_qqq.f 
gel.o: gel.f 
	f77 -c gel.f 
gel_inv.o: gel_inv.f 
	f77 -c gel_inv.f 
prtcl_3d_geo.o: prtcl_3d_geo.f 
	f77 -c prtcl_3d_geo.f 
interp_p.o: interp_p.f
	f77 -c interp_p.f
abc.o: abc.f
	f77 -c abc.f
printel.o: printel.f
	f77 -c printel.f
elm_geom.o: elm_geom.f
	f77 -c elm_geom.f
slp_trgl6.o: slp_trgl6.f
	f77 -c  slp_trgl6.f 
slp_trgl6_sing.o: slp_trgl6_sing.f
	f77 -c  slp_trgl6_sing.f 
slp_trgl3_sing.o: slp_trgl3_sing.f
	f77 -c  slp_trgl3_sing.f 
gauss_leg.o: gauss_leg.f 
	f77 -c gauss_leg.f 
gauss_trgl.o: gauss_trgl.f 
	f77 -c gauss_trgl.f 
#
# clean
# -----
#
clean:
	rm -f core
	rm -f $(OBJ) prtcl_3d_mob
	rm -f prtcl_3d_mob.net prtcl_3d_mob.out
	rm -f matrix_inverse.out
	rm -f particle_elements.out
#
# purge
# ---
#
purge:
	rm -f core 
	rm -f $(OBJ) prtcl_3d_mob
	rm -f prtcl_3d_mob.net prtcl_3d_mob.out
	rm -f matrix_inverse.out
	rm -f particle_elements.out
#
# clobber
# ---
#
clobber:
	rm *
#
# all
# ---
#
all:
	make prtcl_3d_mob
