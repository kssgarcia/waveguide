#
# object files
# ----
#
OBJ1 = trgl6_dr.o trgl6_octa.o trgl6_icos.o
OBJ1A = printel.o printel_stl.o abc.o interp_vn_hs.o elm_geom.o interp_p.o
OBJA = srf_int_3d.o srf_int_3d_interp.o srf_int_3d_ir.o srf_int_3d_ir_interp.o
OBJ2 = cramer_33.o leg_fnc_ortho.o gauss_trgl.o
OBJ3 = vtg_3d.o vtg_3d_interp.o crvm6.o crvm6_interp.o crvten6_interp.o
OBJB = dilt_3d.o dilt_3d_interp.o cmm.o cmm_interp_elm.o int_rule.o
OBJZ = $(OBJ1) $(OBJ1A) $(OBJA) $(OBJ2) $(OBJ3) $(OBJB)
OBJ  = $(OBJZ) 
#
OBJ1hsph = trgl6_hsph_dr.o trgl6_hsph_octa.o
OBJ1Ahsph = printel.o printel_stl.o abc.o interp_vn_hs.o elm_geom.o interp_p.o
OBJ2hsph = gauss_trgl.o crvm6.o crvm6_interp.o crvten6_interp.o
OBJ_hsph = $(OBJ1hsph) $(OBJ1Ahsph) $(OBJ2hsph)
#
# trgl6
# ----
#
trgl6: $(OBJ)
	f77 -o trgl6 $(OBJ)
#
# trgl6_hsph
# --------
#
trgl6_hsph: $(OBJ_hsph)
	f77 -o trgl6_hsph $(OBJ_hsph)
#
# compile:
# -------
#
cramer_33.o: cramer_33.f
	f77 -c cramer_33.f
crvm6.o: crvm6.f
	f77 -c crvm6.f
crvm6_interp.o: crvm6_interp.f
	f77 -c crvm6_interp.f
dilt_3d.o: dilt_3d.f
	f77 -c dilt_3d.f
dilt_3d_interp.o: dilt_3d_interp.f
	f77 -c dilt_3d_interp.f
cmm.o: cmm.f
	f77 -c cmm.f
cmm_interp_elm.o: cmm_interp_elm.f
	f77 -c cmm_interp_elm.f
gauss_trgl.o: gauss_trgl.f
	f77 -c gauss_trgl.f
leg_fnc_ortho.o: leg_fnc_ortho.f
	f77 -c leg_fnc_ortho.f
srf_int_3d.o: srf_int_3d.f
	f77 -c srf_int_3d.f
srf_int_3d_interp.o: srf_int_3d_interp.f
	f77 -c srf_int_3d_interp.f
srf_int_3d_ir.o: srf_int_3d_ir.f
	f77 -c srf_int_3d_ir.f
srf_int_3d_ir_interp.o: srf_int_3d_ir_interp.f
	f77 -c srf_int_3d_ir_interp.f
trgl6_dr.o: trgl6_dr.f
	f77 -c trgl6_dr.f
trgl6_octa.o: trgl6_octa.f
	f77 -c trgl6_octa.f
trgl6_icos.o: trgl6_icos.f
	f77 -c trgl6_icos.f
interp_p.o: interp_p.f
	f77 -c interp_p.f
printel.o: printel.f
	f77 -c printel.f
printel_stl.o: printel_stl.f
	f77 -c printel_stl.f
abc.o: abc.f
	f77 -c abc.f
interp_vn_hs.o: interp_vn_hs.f
	f77 -c interp_vn_hs.f
elm_geom.o: elm_geom.f
	f77 -c elm_geom.f
vtg_3d.o: vtg_3d.f
	f77 -c vtg_3d.f
vtg_3d_interp.o: vtg_3d_interp.f
	f77 -c vtg_3d_interp.f
int_rule.o: int_rule.f
	f77 -c int_rule.f
#
trgl6_hsph_dr.o: trgl6_hsph_dr.f
	f77 -c trgl6_hsph_dr.f
trgl6_hsph_octa.o: trgl6_hsph_octa.f
	f77 -c trgl6_hsph_octa.f
#
# clean:
# -----
#
clean:
	rm -f core 
	rm -f trgl6.net trgl6.stl
	rm -f trgl6_hsph.net trgl6_hsph.stl
	rm -f $(OBJ) trgl6
	rm -f $(OBJ_hsph) trgl6_hsph
#
# purge:
# -----
#
purge:
	rm -f core 
	rm -f trgl6.net trgl6.stl
	rm -f trgl6_hsph.net trgl6_hsph.stl
	rm -f $(OBJ) 
	rm -f $(OBJ_hsph)
#
#
# all:
# ---
#
all:
	make trgl6
	make trgl6_hsph
#
# tar
# ---
#
tar:
	tar cvf trgl6.tar makefile *.f *.m *.dat
