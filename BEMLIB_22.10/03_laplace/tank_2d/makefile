OBJg =
#
OBJ1 = tank_2d.o gauss_leg.o gel.o 
OBJ2 = tank_2d_sdlp.o lgf_2d_www.o 
OBJs = $(OBJ1) $(OBJ2)
OBJ  = $(OBJs) $(OBJg)
#
# libraries
# --------
#
LIBS = 
#
# compiler flags
# --------------
#
FLAGS =
#
# link
# ----
#
tank_2d: $(OBJ)
	f77 -o tank_2d $(OBJ) $(LIBS)
#
# compile engine
# --------------
#
tank_2d.o: tank_2d.f
	f77 -c tank_2d.f
tank_2d_sdlp.o: tank_2d_sdlp.f
	f77 -c tank_2d_sdlp.f
gauss_leg.o: gauss_leg.f
	f77 -c gauss_leg.f
gel.o: gel.f
	f77 -c gel.f
lgf_2d_www.o: lgf_2d_www.f
	f77 -c lgf_2d_www.f
tank_2d_gr.o: tank_2d_gr.f
	f77 -c $(FLAGS) tank_2d_gr.f
#
# all
# -----
#
all:
	make tank_2d
#
# clean
# -----
#
clean:
	rm -f core $(OBJs) xyunit tank_2d tank_2d.xy
	rm -f tank_2d.rst tank_2d
#
# purge
# -----
#
purge:
	rm -f core $(OBJs) xyunit tank_2d tank_2d.xy
	rm -f tank_2d.rst
