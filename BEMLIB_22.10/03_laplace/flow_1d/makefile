OBJ0 = elm_line.o elm_arc.o flow_1d.o gauss_leg.o gel.o
OBJ1 = flow_1d_geo.o flow_1d_sdlp.o flow_1d_vel.o lgf_2d_fs.o
OBJ = $(OBJ0) $(OBJ1)
#
# link
# ----
#
flow_1d: $(OBJ)
	f77 -o flow_1d $(OBJ)
#
# compile
# -------
#
flow_1d.o: flow_1d.f
	f77 -c flow_1d.f
flow_1d_geo.o: flow_1d_geo.f
	f77 -c flow_1d_geo.f
flow_1d_sdlp.o: flow_1d_sdlp.f
	f77 -c flow_1d_sdlp.f
flow_1d_vel.o: flow_1d_vel.f
	f77 -c flow_1d_vel.f
elm_line.o: elm_line.f
	f77 -c elm_line.f
elm_arc.o: elm_arc.f
	f77 -c elm_arc.f
gauss_leg.o: gauss_leg.f
	f77 -c gauss_leg.f
gel.o: gel.f
	f77 -c gel.f
lgf_2d_fs.o: lgf_2d_fs.f
	f77 -c lgf_2d_fs.f
#
# all
# ---
#
all:
	make flow_1d
#
# tar
# ---
#
tar:
	tar cvf flow_1d.tar *.f makefile *.dat
#
# clean
# -----
#
clean:
	rm -f core 
	rm -f $(OBJ) flow_1d
	rm -f flow_1d.elm flow_1d.trc flow_1d.vel flow_1d.profile
#
# purge
# -----
#
purge:
	rm -f core
	rm -f $(OBJ) 
	rm -f flow_1d.elm flow_1d.trc flow_1d.vel flow_1d.profile
