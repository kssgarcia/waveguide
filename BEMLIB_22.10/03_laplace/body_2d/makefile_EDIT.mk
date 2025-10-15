# Fortran compiler
# ----------------
FC = gfortran
FFLAGS =

OBJg =
OBJ1 = elm_line.o elm_arc.o body_2d_EDIT.o gauss_leg.o gel.o
OBJ2 = body_2d_sdlp.o body_2d_vel.o lgf_2d_fs.o lgf_2d_w.o 
OBJ3 = body_2d_geo.o 
OBJs = $(OBJ1) $(OBJ2) $(OBJ3)

OBJ = $(OBJs) $(OBJg)

# libraries
# ---------
LIBS = -lm 

# link
# ----
body_2d_EDIT: $(OBJ)
	$(FC) -o $@ $(OBJ) $(LIBS)

# compile source
# --------------
%.o: %.f
	$(FC) $(FFLAGS) -c $<

# all
# ---
all: body_2d_EDIT

# clean
# -----
clean:
	rm -f core xyunit PLOTDAT
	rm -f $(OBJs) body_2d_EDIT
	rm -f body_2d.str body_2d.out

# purge
# -----
purge:
	rm -f core xyunit PLOTDAT
	rm -f $(OBJs)
	rm -f body_2d.str body_2d.out
