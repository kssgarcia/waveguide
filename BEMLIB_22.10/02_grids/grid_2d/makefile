OBJg =
OBJs = elm_line.o elm_arc.o grid_2d.o
OBJ  = $(OBJs) $(OBJg)
#
# Libraries
# ---------
#
LIBS = 
#
# compiler flags
# --------------
#
FLAGS =
#
#
# link
# ----
#
grid_2d: $(OBJ)
	f77 -o grid_2d $(OBJ) $(LIBS)
#
# compile
# -------
#
grid_2d.o: grid_2d.f
	f77 -c grid_2d.f
elm_line.o: elm_line.f
	f77 -c elm_line.f
elm_arc.o: elm_arc.f
	f77 -c elm_arc.f
#
# clean
# -----
#
clean:
	rm -f core* PLOTDAT
	rm -f $(OBJs) xyunit
	rm -f $(OBJs) grid_2d
#
# purge
# -----
#
purge:
	rm -f core* PLOTDAT
	rm -f $(OBJs) xyunit 
#
# all
# ---
#
all:
	make grid_2d
