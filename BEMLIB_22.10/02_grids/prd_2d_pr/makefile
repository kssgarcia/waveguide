#
# prd_2d_pr
# ---------
#
OBJ1 = prd_2d_pr_dr.o
OBJ2 = prd_2d_pr_arc.o prd_2d_pr_arc_interp.o arc_2d.o 
OBJ3 = prd_2d_pr_splc.o prd_2d_pr_splc_interp.o splc_pr.o splc_pr_geo.o thomas_pr.o
OBJ = $(OBJ1) $(OBJ2) $(OBJ3)
#
#
prd_2d_pr: $(OBJ)
	f77 -o prd_2d_pr $(OBJ)
#
prd_2d_pr_dr.o: prd_2d_pr_dr.f
	f77 -c prd_2d_pr_dr.f
prd_2d_pr_arc.o: prd_2d_pr_arc.f
	f77 -c prd_2d_pr_arc.f
prd_2d_pr_arc_interp.o: prd_2d_pr_arc_interp.f
	f77 -c prd_2d_pr_arc_interp.f
arc_2d.o: arc_2d.f
	f77 -c arc_2d.f
prd_2d_pr_splc.o: prd_2d_pr_splc.f
	f77 -c prd_2d_pr_splc.f
prd_2d_pr_splc_interp.o: prd_2d_pr_splc_interp.f
	f77 -c prd_2d_pr_splc_interp.f
splc_pr.o: splc_pr.f
	f77 -c $(FLAGS) splc_pr.f
splc_pr_geo.o: splc_pr_geo.f
	f77 -c $(FLAGS) splc_pr_geo.f
thomas_pr.o: thomas_pr.f
	f77 -c $(FLAGS) thomas_pr.f
#
# all
# ---
#
all:
	make prd_2d_pr
#
# clean
# -----
#
clean:
	rm -f core 
	rm -f $(OBJ) prd_2d_pr prd_2d_pr.out
