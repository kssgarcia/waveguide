#
# objects
# -------
#
OBJ = prd_3d_pr_dr.o prd_3d_pr.o arc_3d.o cramer_33.o
#
# prd_3d_pr      ! periodic line
# ---------
#
prd_3d_pr: $(OBJ)
	f77 -o prd_3d_pr $(OBJ)
prd_3d_pr.o: prd_3d_pr.f
	f77 -c prd_3d_pr.f
arc_3d.o: arc_3d.f
	f77 -c arc_3d.f
cramer_33.o: cramer_33.f
	f77 -c cramer_33.f
#
# all
# -----
#
all:
	make prd_3d_pr
#
# clean
# -----
#
clean:
	rm -f core
	rm -f $(OBJ) prd_3d_pr prd_3d_pr.net
