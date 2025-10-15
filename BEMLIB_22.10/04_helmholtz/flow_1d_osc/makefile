OBJ1 = flow_1d_osc.o flow_1d_osc_geo.o
OBJ2 = elm_line.o elm_arc.o gauss_leg.o gel.o
OBJ3 = flow_1d_osc_sdlp.o flow_1d_osc_vel.o hgf_2dc_fs.o
OBJ4 = ker_kei_0.o ker_kei_0_large.o
OBJ = $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4)
#
# link
# ----
#
flow_1d_osc: $(OBJ)
	f77 -o flow_1d_osc $(OBJ)
#
# compile
# -------
#
flow_1d_osc.o: flow_1d_osc.f
	f77 -c flow_1d_osc.f
flow_1d_geo.o: flow_1d_geo.f
	f77 -c flow_1d_geo.f
flow_1d_osc_sdlp.o: flow_1d_osc_sdlp.f
	f77 -c flow_1d_osc_sdlp.f
flow_1d_osc_vel.o: flow_1d_osc_vel.f
	f77 -c flow_1d_osc_vel.f
elm_line.o: elm_line.f
	f77 -c elm_line.f
elm_arc.o: elm_arc.f
	f77 -c elm_arc.f
gauss_leg.o: gauss_leg.f
	f77 -c gauss_leg.f
gel.o: gel.f
	f77 -c gel.f
ker_kei_0.o: ker_kei_0.f
	f77 -c ker_kei_0.f
ker_kei_0_large.o: ker_kei_0_large.f
	f77 -c ker_kei_0_large.f
hgf_2dc_fs.o: hgf_2dc_fs.f
	f77 -c hgf_2dc_fs.f
#
# all
# ---
#
all:
	make flow_1d_osc
#
# tar
# ---
#
tar:
	tar cvf flow_1d_osc.tar makefile *.f *dat
#
# clean
#- ----
#
clean:
	rm -f core 
	rm -f $(OBJ) flow_1d_osc
	rm -f flow_1d_osc.trc flow_1d_osc.vel
	rm -f flow_1d_osc.elm
#
# purge
#- ----
#
purge:
	rm -f core
	rm -f $(OBJ) 
	rm -f flow_1d_osc.trc flow_1d_osc.vel
	rm -f flow_1d_osc.elm
