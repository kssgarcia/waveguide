#
# object files
# ------------
#
OBJ1  = ldr_3d.o
OBJ10 = ldlp_3d.o ldlp_3d_interp.o ldlp_3d_integral.o
OBJ11 = lslp_3d.o lslp_3d_interp.o lslp_3d_integral.o lslp_3d_integral_sing.o
OBJ12 = lsdlpp_3d.o lsdlpp_3d_integral.o lsdlpp_3d_interp.o
OBJ2  = lgf_3d_fs.o
OBJ3  = trgl6_octa.o trgl6_icos.o gauss_trgl.o gauss_leg.o
OBJ30 = abc.o interp_p.o elm_geom.o gel.o printel.o
OBJ  = $(OBJ1) $(OBJ10) $(OBJ11) $(OBJ12) $(OBJ2) $(OBJ3) $(OBJ30)
#
# link
# ----
#
ldr_3d: $(OBJ)
	f77 -o ldr_3d $(OBJ)
#
# compile
# -------
#
ldr_3d.o: ldr_3d.f 
	f77 -c ldr_3d.f 
trgl6_octa.o: trgl6_octa.f 
	f77 -c trgl6_octa.f 
trgl6_icos.o: trgl6_icos.f 
	f77 -c trgl6_icos.f 
lgf_3d_fs.o: lgf_3d_fs.f 
	f77 -c lgf_3d_fs.f 
lslp_3d.o: lslp_3d.f
	f77 -c lslp_3d.f 
lslp_3d_interp.o: lslp_3d_interp.f
	f77 -c lslp_3d_interp.f 
lslp_3d_integral.o: lslp_3d_integral.f
	f77 -c lslp_3d_integral.f 
lslp_3d_integral_sing.o: lslp_3d_integral_sing.f
	f77 -c lslp_3d_integral_sing.f 
ldlp_3d.o: ldlp_3d.f
	f77 -c ldlp_3d.f 
ldlp_3d_interp.o: ldlp_3d_interp.f
	f77 -c ldlp_3d_interp.f 
ldlp_3d_integral.o: ldlp_3d_integral.f
	f77 -c ldlp_3d_integral.f 
lsdlpp_3d.o: lsdlpp_3d.f 
	f77 -c lsdlpp_3d.f 
lsdlpp_3d_interp.o: lsdlpp_3d_interp.f 
	f77 -c lsdlpp_3d_interp.f 
lsdlpp_3d_integral.o: lsdlpp_3d_integral.f 
	f77 -c lsdlpp_3d_integral.f 
gel.o: gel.f
	f77 -c gel.f 
abc.o: abc.f
	f77 -c abc.f
elm_geom.o: elm_geom.f 
	f77 -c elm_geom.f 
interp_p.o: interp_p.f
	f77 -c interp_p.f
printel.o: printel.f
	f77 -c printel.f
gauss_leg.o: gauss_leg.f 
	f77 -c gauss_leg.f 
gauss_trgl.o: gauss_trgl.f 
	f77 -c gauss_trgl.f 
#
# all
# ---
#
all:
	make ldr_3d
#
# tar
# ---
#
tar:
	tar cvf ldr_3d.tar makefile *.f *.dat
#
# clean
# -----
#
clean:
	rm -f core
	rm -f $(OBJ) ldr_3d
	rm -f ldr_3d.out ldr_3d.net
#
# purge
# -----
#
purge:
	rm -f core
	rm -f $(OBJ)
	rm -f ldr_3d.out ldr_3d.net
