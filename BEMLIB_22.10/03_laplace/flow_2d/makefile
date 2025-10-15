OBJg =
#
OBJ1 = elm_line.o elm_arc.o flow_2d.o gauss_leg.o gel.o
OBJ2 = flow_2d_sdlp.o flow_2d_vel.o lgf_2d_fs.o 
OBJ3 = flow_2d_geo.o
OBJs = $(OBJ1) $(OBJ2) $(OBJ3)
OBJ  = $(OBJs) $(OBJg)
#
# Libraries
# ---------
#
LIBS = 
#
# compiler flags
# ---------
#
FLAGS =
#
# link
# ----
#
flow_2d: $(OBJ)
	f77 -o flow_2d $(OBJ) $(LIBS)
#
# compile source
# --------------
#
flow_2d.o: flow_2d.f
	f77 -c flow_2d.f
flow_2d_geo.o: flow_2d_geo.f
	f77 -c flow_2d_geo.f
flow_2d_sdlp.o: flow_2d_sdlp.f
	f77 -c flow_2d_sdlp.f
flow_2d_vel.o: flow_2d_vel.f
	f77 -c flow_2d_vel.f
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
# -----
#
all:
	make flow_2d
#
# clean
# -----
#
clean:
	rm -f core xyunit PLOTDAT
	rm -f $(OBJs) flow_2d
	rm -f flow_2d.prof flow_2d.str flow_2d.out
#
# purge
# -----
#
purge:
	rm -f core xyunit PLOTDAT
	rm -f $(OBJs)
	rm -f flow_2d.prof flow_2d.str flow_2d.out
