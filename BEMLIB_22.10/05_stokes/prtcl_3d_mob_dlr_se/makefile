OBJ1 = prtcl_3d_mob_dlr_se.o prtcl_3d_mob_dlr_se_slv.o
OBJ2 = sgf_3d_fs.o sgf_3d_w.o
OBJ3 = trgl6_octa.o trgl6_icos.o
OBJ9 = gauss_trgl.o
OBJ4 = jacobi.o gel_inv.o
OBJ8 = sdlp_3d.o sdlp_3d_integral.o sdlp_3d_interp.o 
OBJ5 = vander.o lobatto.o trgl_rule.o
OBJ30 = elm_geom.o abc.o interp_p.o printel.o
OBJ  = $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ8) $(OBJ9) $(OBJ30) 
#
# link
# ----
#
prtcl_3d_mob_dlr_se: $(OBJ)
	f77 -o prtcl_3d_mob_dlr_se $(OBJ)
#
# compile
# -------
#
prtcl_3d_mob_dlr_se.o: prtcl_3d_mob_dlr_se.f 
	f77 -c prtcl_3d_mob_dlr_se.f 
prtcl_3d_mob_dlr_se_slv.o: prtcl_3d_mob_dlr_se_slv.f 
	f77 -c prtcl_3d_mob_dlr_se_slv.f 
printel.o: printel.f
	f77 -c printel.f
abc.o: abc.f
	f77 -c abc.f
interp_p.o: interp_p.f
	f77 -c interp_p.f
elm_geom.o: elm_geom.f
	f77 -c elm_geom.f
trgl6_octa.o: trgl6_octa.f 
	f77 -c trgl6_octa.f 
trgl6_icos.o: trgl6_icos.f 
	f77 -c trgl6_icos.f 
trgl_octa.o: trgl_octa.f 
	f77 -c trgl_octa.f 
sgf_3d_fs.o: sgf_3d_fs.f 
	f77 -c sgf_3d_fs.f 
sgf_3d_w.o: sgf_3d_w.f 
	f77 -c sgf_3d_w.f 
sdlp_3d.o: sdlp_3d.f
	f77 -c sdlp_3d.f 
sdlp_3d_interp.o: sdlp_3d_interp.f
	f77 -c sdlp_3d_interp.f 
sdlp_3d_integral.o: sdlp_3d_integral.f
	f77 -c sdlp_3d_integral.f 
gauss_trgl.o: gauss_trgl.f
	f77 -c gauss_trgl.f 
vander.o: vander.f
	f77 -c vander.f 
lobatto.o: lobatto.f
	f77 -c lobatto.f 
jacobi.o: jacobi.f
	f77 -c jacobi.f 
gel_inv.o: gel_inv.f
	f77 -c gel_inv.f 
trgl_rule.o: trgl_rule.f
	f77 -c trgl_rule.f 
#
# all
# ---
#
all:
	make prtcl_3d_mob_dlr_se
#
# clean
# -----
#
clean:
	rm -f core 
	rm -f $(OBJ) prtcl_3d_mob_dlr_se
	rm -f prtcl_3d_mob_dlr_se.out prtcl_3d_mob_dlr_se.net
#
# purge
# -----
#
purge:
	rm -f core
	rm -f $(OBJ)
	rm -f prtcl_3d_mob_dlr_se.out prtcl_3d_mob_dlr_se.net
