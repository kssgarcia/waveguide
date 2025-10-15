#
# Object files
#
OBJ1 = ldr_3d_2p.o trgl6_sqr.o gauss_trgl.o gauss_leg.o
OBJ10 = ldlp_3d_2p.o ldlp_3d_interp.o ldlp_3d_2p_integral.o
OBJ12 = ldlpp_3d_2p.o ldlpp_3d_2p_interp.o ldlpp_3d_2p_integral.o
OBJ2 = lslp3_3d_2p.o
OBJ3 = lgf_3d_2p.o
OBJ30 = abc.o interp_p.o elm_geom.o printel.o
OBJ4 = sgrad_3d_2p.o sgrad_3d_2p_interp.o cramer_33.o
OBJ41 = vs_3d_2p_circ_vort.o vs_3d_2p_circ_vort_interp.o vs_3d_2p_curl.o vs_3d_2p_curl_interp.o 
OBJ5 = srf_int_3d_2p.o
OBJ  = $(OBJ1) $(OBJ10) $(OBJ12) $(OBJ2) $(OBJ3) $(OBJ30) $(OBJ4) $(OBJ41) $(OBJ5)
#
# link
# ----
#
ldr_3d_2p: $(OBJ)
	f77 -o ldr_3d_2p $(OBJ)
#
# compile
# --------
#
ldr_3d_2p.o: ldr_3d_2p.f 
	f77 -c ldr_3d_2p.f 
ldr_3d_2p_geo.o: ldr_3d_2p_geo.f 
	f77 -c ldr_3d_2p_geo.f 
abc.o: abc.f
	f77 -c abc.f
interp_p.o: interp_p.f
	f77 -c interp_p.f
elm_geom.o: elm_geom.f
	f77 -c elm_geom.f
printel.o: printel.f
	f77 -c printel.f
trgl_sqr.o: trgl_sqr.f 
	f77 -c trgl_sqr.f 
lgf_3d_2p.o: lgf_3d_2p.f 
	f77 -c lgf_3d_2p.f 
ldlp_3d_2p.o: ldlp_3d_2p.f
	f77 -c ldlp_3d_2p.f 
ldlp_3d_interp.o: ldlp_3d_interp.f
	f77 -c ldlp_3d_interp.f 
ldlp_3d_2p_integral.o: ldlp_3d_2p_integral.f
	f77 -c ldlp_3d_2p_integral.f 
ldlpp_3d_2p.o: ldlpp_3d_2p.f
	f77 -c ldlpp_3d_2p.f 
ldlpp_3d_2p_interp.o: ldlpp_3d_2p_interp.f
	f77 -c ldlpp_3d_2p_interp.f 
ldlpp_3d_2p_integral.o: ldlpp_3d_2p_integral.f
	f77 -c ldlpp_3d_2p_integral.f 
gauss_trgl.o: gauss_trgl.f
	f77 -c gauss_trgl.f 
gauss_leg.o: gauss_leg.f
	f77 -c gauss_leg.f 
lslp3_3d_2p.o: lslp3_3d_2p.f
	f77 -c lslp3_3d_2p.f 
sgrad_3d_2p.o: sgrad_3d_2p.f
	f77 -c sgrad_3d_2p.f 
sgrad_3d_2p_interp.o: sgrad_3d_2p_interp.f
	f77 -c sgrad_3d_2p_interp.f 
vs_3d_2p_circ_vort.o: vs_3d_2p_circ_vort.f
	f77 -c vs_3d_2p_circ_vort.f 
vs_3d_2p_circ_vort_interp.o: vs_3d_2p_circ_vort_interp.f
	f77 -c vs_3d_2p_circ_vort_interp.f 
vs_3d_2p_curl.o: vs_3d_2p_curl.f
	f77 -c vs_3d_2p_curl.f 
vs_3d_2p_curl_interp.o: vs_3d_2p_curl_interp.f
	f77 -c vs_3d_2p_curl_interp.f 
cramer_33.o: cramer_33.f
	f77 -c cramer_33.f 
srf_int_3d_2p.o: srf_int_3d_2p.f
	f77 -c srf_int_3d_2p.f
#
# all
# ---
#
all:
	make ldr_3d_2p
#
#
# clean
# -----
#
clean:
	rm -f core 
	rm -f $(OBJ) ldr_3d_2p
	rm -f ldr_3d_2p ldr_3d_2p.sln ldr_3d_2p.grid ldr_3d_2p.net ldr_3d_2p.out
