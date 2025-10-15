OBJg =
OBJ1 = elm_line.o elm_arc.o body_ax.o gauss_leg.o gel.o
OBJ2 = body_ax_sdlp.o body_ax_vel.o lgf_ax_fs.o ell_int.o
OBJ3 = lvr_fs.o body_ax_geo.o
OBJs = $(OBJ1) $(OBJ2) $(OBJ3)
OBJ  = $(OBJs) $(OBJg)
#
# Libraries
# ---------
#
LIBS = -lm
#
# compiler flags
# ---------
#
FLAGS =
#
# link
# -----
#
body_ax: $(OBJ)
	f77 -o body_ax $(OBJ) $(LIBS)
#
# compile
# -------
#
body_ax.o: body_ax.f
	f77 -c body_ax.f
body_ax_sdlp.o: body_ax_sdlp.f
	f77 -c body_ax_sdlp.f
body_ax_geo.o: body_ax_geo.f
	f77 -c body_ax_geo.f
body_ax_vel.o: body_ax_vel.f
	f77 -c body_ax_vel.f
lvr_fs.o: lvr_fs.f
	f77 -c lvr_fs.f
elm_line.o: elm_line.f
	f77 -c elm_line.f
elm_arc.o: elm_arc.f
	f77 -c elm_arc.f
gauss_leg.o: gauss_leg.f
	f77 -c gauss_leg.f
gel.o: gel.f
	f77 -c gel.f
lgf_ax_fs.o: lgf_ax_fs.f
	f77 -c lgf_ax_fs.f
ell_int.o: ell_int.f
	f77 -c ell_int.f
#
# all
# -----
#
all:
	make body_ax
#
# clean
# -----
#
clean:
	rm -f core PLOTDAT AAA-xmgr
	rm -f $(OBJs) body_ax
	rm -f body_ax.out body_ax.str xyunit
#
# purge
# -----
#
purge:
	rm -f core PLOTDAT AAA-xmgr
	rm -f $(OBJs) 
	rm -f body_ax.out body_ax.str xyunit
