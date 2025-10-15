OBJ1 = ldr_3d_int.o
OBJ10 = ldlp_3d.o ldlp_3d_interp.o ldlp_3d_integral.o
OBJ11A = lslp3_3d.o lslp3_3d_interp.o lslp3_3d_integrate.o
OBJ11B = lslp3_3d_integrate_sing.o
OBJ12 = ldlpp_3d.o ldlpp_3d_interp.o ldlpp_3d_integral.o
OBJ2 = lgf_3d_fs.o lgf_3d_w.o
OBJ3 = trgl6_octa.o trgl6_icos.o gauss_trgl.o gauss_leg.o
OBJ33 = srf_int_3d.o srf_int_3d_interp.o
OBJ30 = abc.o interp_p.o elm_geom.o printel.o 
OBJ4 = sgrad_3d.o sgrad_3d_interp.o cramer_33.o
OBJ41 = vs_3d_circ_vort.o vs_3d_circ_vort_interp.o vs_3d_curl.o vs_3d_curl_interp.o
OBJA  = $(OBJ1) $(OBJ10) $(OBJ11A) $(OBJ11B) $(OBJ12) $(OBJ2) $(OBJ3)
OBJB  = $(OBJ30) $(OBJ4) $(OBJ41) $(OBJ33)
OBJ   = $(OBJA) $(OBJB)
#
# link
# ----
#
ldr_3d_int: $(OBJ)
	f77 -o ldr_3d_int $(OBJ)
#
# compile
# -------
#
ldr_3d_int.o: ldr_3d_int.f 
	f77 -c ldr_3d_int.f 
ldr_3d_geo.o: ldr_3d_geo.f 
	f77 -c ldr_3d_geo.f 
trgl6_octa.o: trgl6_octa.f 
	f77 -c trgl6_octa.f 
trgl6_icos.o: trgl6_icos.f 
	f77 -c trgl6_icos.f 
lgf_3d_fs.o: lgf_3d_fs.f 
	f77 -c lgf_3d_fs.f 
lgf_3d_w.o: lgf_3d_w.f 
	f77 -c lgf_3d_w.f 
ldlp_3d.o: ldlp_3d.f
	f77 -c ldlp_3d.f 
ldlp_3d_interp.o: ldlp_3d_interp.f
	f77 -c ldlp_3d_interp.f 
ldlp_3d_integral.o: ldlp_3d_integral.f
	f77 -c ldlp_3d_integral.f 
ldlpp_3d.o: ldlp_3d.f
	f77 -c ldlpp_3d.f 
ldlpp_3d_interp.o: ldlpp_3d_interp.f
	f77 -c ldlpp_3d_interp.f 
ldlpp_3d_integral.o: ldlpp_3d_integral.f
	f77 -c ldlpp_3d_integral.f 
gauss_trgl.o: gauss_trgl.f
	f77 -c gauss_trgl.f 
gauss_leg.o: gauss_leg.f
	f77 -c gauss_leg.f 
lslp3_3d.o: lslp3_3d.f
	f77 -c lslp3_3d.f 
lslp3_3d_interp.o: lslp3_3d_interp.f
	f77 -c lslp3_3d_interp.f 
lslp3_3d_integrate.o: lslp3_3d_integrate.f
	f77 -c lslp3_3d_integrate.f 
lslp3_3d_integrate_sing.o: lslp3_3d_integrate_sing.f
	f77 -c lslp3_3d_integrate_sing.f 
abc.o: abc.f
	f77 -c abc.f
interp_p.o: interp_p.f
	f77 -c interp_p.f
elm_geom.o: elm_geom.f
	f77 -c elm_geom.f
printel.o: printel.f
	f77 -c printel.f
sgrad_3d.o: sgrad_3d.f
	f77 -c sgrad_3d.f 
sgrad_3d_interp.o: sgrad_3d_interp.f
	f77 -c sgrad_3d_interp.f 
vs_3d_circ_vort.o: vs_3d_circ_vort.f
	f77 -c vs_3d_circ_vort.f 
vs_3d_circ_vort_interp.o: vs_3d_circ_vort_interp.f
	f77 -c vs_3d_circ_vort_interp.f 
vs_3d_curl.o: vs_3d_curl.f
	f77 -c vs_3d_curl.f 
vs_3d_curl_interp.o: vs_3d_curl_interp.f
	f77 -c vs_3d_curl_interp.f 
cramer_33.o: cramer_33.f
	f77 -c cramer_33.f 
srf_int_3d.o: srf_int_3d.f
	f77 -c srf_int_3d.f
srf_int_3d_interp.o: srf_int_3d_interp.f
	f77 -c srf_int_3d_interp.f
#
# clean
# -----
#
clean:
	rm -f core 
	rm -f $(OBJ) ldr_3d_int
	rm -f ldr_3d_int.out ldr_3d_int.net ldr_3d_int_1.net
#
# purge
# -----
#
purge:
	rm -f core 
	rm -f $(OBJ)
	rm -f ldr_3d_int.out ldr_3d_int.net ldr_3d_int_1.net
# all
# ---
#
all:
	make ldr_3d_int
#
# mpage
# -----
#
mult_page:
	mpage -4 -c -bLetter makefile ldr_3d_int.f geo_3d.f ldlp_3d.f ldlpp_3d.f lslp3_3d.f lgf_3d_fs.f vs_3d_circ_vort.f vs_3d_curl.f cramer_33.f trgl_octa.f gauss_trgl.f gauss_leg.f srf_int_3d.f sgrad_3d.f ldr_3d_int.dat > tmp
