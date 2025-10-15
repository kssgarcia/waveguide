#
# objects:
# -------
#
OBJg =
#
OBJs = drop_ax.o drop_ax_ode.o
OBJ  = $(OBJs) $(OBJg)
#
# libraries:
# ---------
#
LIBS =
#
# compiler flags:
# --------------
#
FLAGS =
#
# link:
# ----
#
drop_ax: $(OBJ)
	f77 -o drop_ax $(OBJ) $(LIBS)
#
# compile:
# -------
#
drop_ax.o: drop_ax.f
	f77 -c $(FLAGS) drop_ax.f
drop_ax_ode.o: drop_ax_ode.f
	f77 -c $(FLAGS) drop_ax_ode.f
#
# ode
# -----
#
ode: ode.o ode_vel.o
	f77 -o ode ode.o ode_vel.o
ode.o: ode.f
	f77 -c ode.f
ode_vel.o: ode_vel.f
	f77 -c ode_vel.f
#
# all:
# ---
#
all:
	make drop_ax
	make ode
#
# clean:
# -----
#
clean:
	rm -f core $(OBJs) xyunit drop_ax.out
	rm -f drop_ax
	rm -f ode.o ode_vel.o ode ode.out
