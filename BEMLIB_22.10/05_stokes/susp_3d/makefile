#
# graphics:
# ---------
#
OBJg =
#
#
# object files
# ------------
#
OBJ1 = susp_3d.o susp_3d_col.o susp_3d_slv.o
OBJ2 = susp_3d_dlp.o susp_3d_dlp_integral.o susp_3d_dlp_interp.o
OBJ2A = susp_3d_dlp_mutual.o susp_3d_dlp_integral_mutual.o
OBJ3 = sgf_3d_fs.o sgf_3d_w.o
OBJ31 = sgf_3d_3p.o sgf_3d_3p_qqq.o sgf_3d_3p_vvv.o sgf_3d_3p_ewald.o
OBJ30 = elm_geom.o abc.o interp_p.o printel.o
OBJ4 = gauss_trgl.o vander.o lobatto.o trgl_rule.o
OBJ41 = trgl6_octa.o trgl6_icos.o
OBJ5 = jacobi.o gel_inv.o cramer_33.o
OBJs = $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ30) $(OBJ31) $(OBJ41) $(OBJ2A)
#
OBJ = $(OBJs) $(OBJg)
#
# libraries
# ---------
#
LIBS = 
#
# link
# ----
#
susp_3d: $(OBJ)
	f77 -o susp_3d $(OBJ) $(LIBS)
prytania: $(OBJs)
	f77 -o prytania $(OBJs)
#
# compile
# -------
#
cramer_33.o: cramer_33.f 
	f77 -c cramer_33.f 
susp_3d.o: susp_3d.f 
	f77 -c susp_3d.f 
susp_3d_slv.o: susp_3d_slv.f 
	f77 -c susp_3d_slv.f 
printel.o: printel.f
	f77 -c printel.f
abc.o: abc.f
	f77 -c abc.f
interp_p.o: interp_p.f
	f77 -c interp_p.f
elm_geom.o: elm_geom.f
	f77 -c elm_geom.f
susp_3d_dlp.o: susp_3d_dlp.f
	f77 -c susp_3d_dlp.f
susp_3d_dlp_integral.o: susp_3d_dlp_integral.f
	f77 -c susp_3d_dlp_integral.f
susp_3d_dlp_interp.o: susp_3d_dlp_interp.f
	f77 -c susp_3d_dlp_interp.f
susp_3d_dlp_mutual.o: susp_3d_dlp_mutual.f
	f77 -c susp_3d_dlp_mutual.f
susp_3d_dlp_integral_mutual.o: susp_3d_dlp_integral_mutual.f
	f77 -c susp_3d_dlp_integral_mutual.f
susp_3d_col.o: susp_3d_col.f
	f77 -c susp_3d_col.f
trgl6_octa.o: trgl6_octa.f 
	f77 -c trgl6_octa.f 
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
sgf_3d_3p_vvv.o: sgf_3d_3p_vvv.f 
	f77 -c sgf_3d_3p_vvv.f 
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
	make susp_3d
#
# tar
# ---
#
tar:
	tar cvf susp_3d.tar makefile *.f *.dat *.m
#
# clean
# -----
#
clean:
	rm -f core
	rm -f $(OBJ) susp_3d
	rm -f susp_3d.rhe susp_3d.trj susp_3d.orb susp_3d.adh susp_3d.vel
	rm -f trgl_3d.net
