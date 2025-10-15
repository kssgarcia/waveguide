#
# prd_ax
# ------
#
OBJ = prd_ax_dr.o prd_ax.o arc_ax.o sm_ax.o
#
prd_ax: $(OBJ)
	f77 -o prd_ax $(OBJ)
#
#
arc_ax: arc_ax.f
	f77 -o arc_ax .f
prd_ax_dr.o: prd_ax_dr.f
	f77 -c prd_ax_dr.f
prd_ax.o: prd_ax.f
	f77 -c prd_ax.f
sm_ax.o: sm_ax.f
	f77 -c sm_ax.f
#
# all
# ---
#
all:
	make prd_ax
#
# clean
# -----
#
clean:
	rm -f core PLOTDAT
	rm -f $(OBJ) prd_ax.out prd_ax 
#
# purge
# -----
#
purge:
	rm -f core PLOTDAT
	rm -f $(OBJ) prd_ax.out
