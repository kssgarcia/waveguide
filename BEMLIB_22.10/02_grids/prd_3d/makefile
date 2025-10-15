#
# objects
# -------
#
OBJ = prd_3d_dr.o prd_3d.o arc_3d.o cramer_33.o
#
# prd_3d      
# ------
#
prd_3d: $(OBJ)
	f77 -o prd_3d $(OBJ)
prd_3d.o: prd_3d.f
	f77 -c prd_3d.f
arc_3d.o: arc_3d.f
	f77 -c arc_3d.f
cramer_33.o: cramer_33.f
	f77 -c cramer_33.f
#
# all
# ---
#
all:
	make prd_3d
#
# clean
# -----
#
clean:
	rm -f core 
	rm -f $(OBJ) prd_3d prd_3d.net
