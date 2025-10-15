#
# prtcl_sw
# --------
#
OBJ0 = elm_line.o elm_arc.o prtcl_sw.o gauss_leg.o gel.o
OBJ1 = prtcl_sw_geo.o prtcl_sw_slp.o  ell_int.o
OBJ  = $(OBJ0) $(OBJ1)
#
# link
# ----
#
prtcl_sw: $(OBJ)
	f77 -o prtcl_sw $(OBJ)
#
# compile
#
prtcl_sw.o: prtcl_sw.f
	f77 -c prtcl_sw.f
prtcl_sw_slp.o: prtcl_sw_slp.f
	f77 -c prtcl_sw_slp.f
prtcl_sw_geo.o: prtcl_sw_geo.f
	f77 -c prtcl_sw_geo.f
elm_line.o: elm_line.f
	f77 -c elm_line.f
elm_arc.o: elm_arc.f
	f77 -c elm_arc.f
gauss_leg.o: gauss_leg.f
	f77 -c gauss_leg.f
gel.o: gel.f
	f77 -c gel.f
ell_int.o: ell_int.f
	f77 -c ell_int.f
#
# all
# ----
#
all:
	make prtcl_sw
#
# clean
# ----
#
clean:
	rm -f core PLOTDAT
	rm -f $(OBJ) prtcl_sw
	rm -f $(OBJ) prtcl_sw.xy prtcl_sw.trc
#
# purge
# ----
#
purge:
	rm -f core PLOTDAT
	rm -f $(OBJ)
	rm -f $(OBJ) prtcl_sw.xy prtcl_sw.trc
#
# tar
# ----
#
tar:
	tar cvf prtcl_sw.tar makefile *.f *.dat
