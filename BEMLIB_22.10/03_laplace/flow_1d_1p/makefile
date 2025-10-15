OBJ0 = elm_line.o elm_arc.o flow_1d_1p.o flow_1d_1p_geo.o
OBJ1 = gauss_leg.o gel.o
OBJ2 = flow_1d_1p_slp.o flow_1d_1p_vel.o lgf_2d_1p.o
OBJ = $(OBJ0) $(OBJ1) $(OBJ2)
#
# link
# -----
#
flow_1d_1p: $(OBJ)
	f77 -o flow_1d_1p $(OBJ)
#
# compile
# -----
#
flow_1d_1p.o: flow_1d_1p.f
	f77 -c flow_1d_1p.f
flow_1d_1p_geo.o: flow_1d_1p_geo.f
	f77 -c flow_1d_1p_geo.f
flow_1d_1p_slp.o: flow_1d_1p_slp.f
	f77 -c flow_1d_1p_slp.f
flow_1d_1p_vel.o: flow_1d_1p_vel.f
	f77 -c flow_1d_1p_vel.f
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
# clean
# -----
#
clean:
	rm -f core PLOTDAT AAA-xmgr
	rm -f $(OBJ) flow_1d_1p
	rm -f flow_1d_1p.elm flow_1d_1p.trc
	rm -f flow_1d_1p.vel
#
# purge
# -----
#
purge:
	rm -f core PLOTDAT AAA-xmgr
	rm -f $(OBJ)
	rm -f flow_1d_1p.elm flow_1d_1p.trc
	rm -f flow_1d_1p.vel
#
# all
# -----
#
all:
	make flow_1d_1p
#
# tar
# ---
#
tar:
	tar cvf flow_1d_1p.tar *
