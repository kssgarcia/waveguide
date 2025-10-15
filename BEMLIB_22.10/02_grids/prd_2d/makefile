OBJ = prd_2d_dr.o prd_2d.o prd_2d_interp.o arc_2d.o 
#
# prd_2d
# ------
#
prd_2d: $(OBJ)
	f77 -o prd_2d $(OBJ)
#
# compile:
# ------
# 
prd_2d_dr.o: prd_2d_dr.f
	f77 -c prd_2d_dr.f
prd_2d.o: prd_2d.f
	f77 -c prd_2d.f
prd_2d_interp.o: prd_2d_interp.f
	f77 -c prd_2d_interp.f
arc_2d.o: arc_2d.f
	f77 -c arc_2d.f
#
# all
# ---
#
all:
	make prd_2d
#
# clean
# -----
#
clean:
	rm -f core 
	rm -f $(OBJ) prd_2d.out prd_2d
#
# purge
# -----
#
purge:
	rm -f core
	rm -f $(OBJ) prd_2d.out
