#
# object files
# ------------
#
OBJ1 = elm_line.o elm_arc.o flow_3x.o gauss_leg.o
OBJ2 = flow_3x_slp.o sgf_3x.o gel.o flow_3x_vel.o davis.o
OBJ3 = ell_int.o
OBJ  = $(OBJ1) $(OBJ2) $(OBJ3)
#
# link
# ----
#
flow_3x: $(OBJ)
	f77 -o flow_3x $(OBJ)
#
# compile
# -------
#
flow_3x.o: flow_3x.f
	f77 -c flow_3x.f
flow_3x_slp.o: flow_3x_slp.f
	f77 -c flow_3x_slp.f
flow_3x_vel.o: flow_3x_vel.f
	f77 -c flow_3x_vel.f
elm_line.o: elm_line.f
	f77 -c elm_line.f
elm_arc.o: elm_arc.f
	f77 -c elm_arc.f
gauss_leg.o: gauss_leg.f
	f77 -c gauss_leg.f
gel.o: gel.f
	f77 -c gel.f
ell_int: ell_int.f
	f77 -c ell_int.f
sgf_3x.o: sgf_3x.f
	f77 -c sgf_3x.f
davis.o: davis.f
	f77 -c davis.f
#
# clean
# -----
#
clean:
	rm -f $(OBJ) flow_3x
	rm -f flow_3x.out1
	rm -f flow_3x.out2
	rm -f flow_3x.out3
#
# purge
# -----
#
purge:
	rm -f $(OBJ) 
	rm -f flow_3x.out1
	rm -f flow_3x.out2
	rm -f flow_3x.out3
#
# all
# -----
#
all:
	make flow_3x
