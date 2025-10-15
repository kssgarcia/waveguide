OBJg =
#
OBJ1 = elm_line.o elm_arc.o body_2d.o gauss_leg.o gel.o
OBJ2 = body_2d_sdlp.o body_2d_vel.o lgf_2d_fs.o lgf_2d_w.o 
OBJ3 = body_2d_geo.o 
OBJs = $(OBJ1) $(OBJ2) $(OBJ3)
#
OBJ = $(OBJs) $(OBJg)
#
# libraries
# ---------
#
LIBS = -lm 
#
# compiler flags
# --------------
#
FLAGS =
#
# link
# ----
#
body_2d: $(OBJ)
	f77 -o body_2d $(OBJ) $(LIBS)
#
# compile source
# --------------
#
body_2d.o: body_2d.f
	f77 -c body_2d.f
body_2d_geo.o: body_2d_geo.f
	f77 -c body_2d_geo.f
body_2d_sdlp.o: body_2d_sdlp.f
	f77 -c body_2d_sdlp.f
body_2d_vel.o: body_2d_vel.f
	f77 -c body_2d_vel.f
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
lgf_2d_w.o: lgf_2d_w.f
	f77 -c lgf_2d_w.f
body_2d_gr.o: body_2d_gr.f
	f77 -c body_2d_gr.f
#
# all
# ---
all:
	make body_2d
#
# clean
# -----
#
clean:
	rm -f core xyunit PLOTDAT
	rm -f $(OBJs) body_2d
	rm -f body_2d.str body_2d.out
#
# purge
# -----
#
purge:
	rm -f core xyunit PLOTDAT
	rm -f $(OBJs)
	rm -f body_2d.str body_2d.out
