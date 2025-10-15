#
# Objects
# -------
#
OBJ0 = verbal.o
OBJ1 = trgl6_octa.o trgl6_octa_ss.o trgl6_icos.o gauss_trgl.o gauss_leg.o
OBJ2 = sgf_3d_fs.o sgf_3d_w.o sgf_3d_2p_w.o
OBJ2A = sgf_3d_3p.o sgf_3d_3p_ewald.o sgf_3d_3p_qqq.o
OBJ2B = sgf_3d_2p.o sgf_3d_2p_ewald.o sgf_3d_2p_aux.o 
OBJ3 = prtcl_3d.o prtcl_3d_prec.o
OBJ30 = elm_geom.o abc.o interp_p.o printel.o
OBJ33 = slp_trgl6.o slp_trgl6_sing.o slp_trgl3_sing.o
OBJ4 = gel_inv.o gel_mrhs.o
OBJ5 = gel_power.o cond_num.o power.o
OBJs = $(OBJ0) $(OBJ1) $(OBJ2) $(OBJ2A) $(OBJ2B) $(OBJ3) $(OBJ30) $(OBJ33) $(OBJ4) $(OBJ5) 
#
# link
# ----
#
prtcl_3d: $(OBJs)
	f77 -o prtcl_3d $(OBJs)
#
# compile
# ------
#
cond_num.o: cond_num.f 
	f77 -c cond_num.f 
gauss_leg.o: gauss_leg.f 
	f77 -c gauss_leg.f 
gauss_trgl.o: gauss_trgl.f 
	f77 -c gauss_trgl.f 
trgl6_octa.o: trgl6_octa.f 
	f77 -c trgl6_octa.f 
trgl6_octa_ss.o: trgl6_octa_ss.f 
	f77 -c trgl6_octa_ss.f 
trgl6_icos.o: trgl6_icos.f 
	f77 -c trgl6_icos.f 
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
sgf_3d_2p.o: sgf_3d_2p.f 
	f77 -c sgf_3d_2p.f 
sgf_3d_2p_w.o: sgf_3d_2p_w.f 
	f77 -c sgf_3d_2p_w.f 
sgf_3d_2p_aux.o: sgf_3d_2p_aux.f 
	f77 -c sgf_3d_2p_aux.f 
sgf_3d_2p_ewald.o: sgf_3d_2p_ewald.f 
	f77 -c sgf_3d_2p_ewald.f 
gel_inv.o: gel_inv.f 
	f77 -c gel_inv.f 
gel_mrhs.o: gel_mrhs.f 
	f77 -c gel_mrhs.f 
gel_power.o: gel_power.f 
	f77 -c gel_power.f 
prtcl_3d.o: prtcl_3d.f 
	f77 -c prtcl_3d.f 
printel.o: printel.f
	f77 -c printel.f
abc.o: abc.f
	f77 -c abc.f
interp_p.o: interp_p.f
	f77 -c interp_p.f
elm_geom.o: elm_geom.f
	f77 -c elm_geom.f
slp_trgl6.o: slp_trgl6.f 
	f77 -c  slp_trgl6.f 
slp_trgl6_sing.o: slp_trgl6_sing.f 
	f77 -c  slp_trgl6_sing.f 
slp_trgl3_sing.o: slp_trgl3_sing.f 
	f77 -c  slp_trgl3_sing.f 
prtcl_3d_prec.o: prtcl_3d_prec.f 
	f77 -c prtcl_3d_prec.f 
power.o: power.f 
	f77 -c power.f 
verbal.o: verbal.f 
	f77 -c verbal.f 
#
# spheroid
# --------
#
spheroid: spheroid.f
	f77 -o spheroid spheroid.f
#
# all
# ---
#
all:
	make prtcl_3d
	make spheroid
#
# clean
# -----
#
clean:
	rm -f core PLOTDAT
	rm -f $(OBJs) prtcl_3d
	rm -f prtcl_3d.out particle_points.out
	rm -f particle_elements.out
	rm -f prtcl_3d.net
	rm -f matrix_inverse.out
	rm -f spheroid.o spheroid
#
# purge
# -----
#
purge:
	rm -f core PLOTDAT
	rm -f prtcl_3d.out
	rm -f $(OBJs)
	rm -f prtcl_3d.out particle_points.out
	rm -f particle_elements.out
	rm -f prtcl_3d.net
	rm -f matrix_inverse.out
	rm -f spheroid.o
