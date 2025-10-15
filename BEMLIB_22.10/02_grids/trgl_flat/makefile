#
# object files
# ------------
#
OBJ0 = trgl6.o trgl6_sqr.o trgl6_sc.o trgl6_ss.o trgl6_disk.o
OBJ1 = trgl6_sqr_dr.o trgl6_sqr.o
OBJ2 = trgl6_sqr_geo.o gauss_trgl.o
OBJ3 = crvm_3d_2p.o srf_int_3d_2p.o
OBJA  = $(OBJ1) $(OBJ2) $(OBJ3)
OBJB  = $(OBJ0) $(OBJ2) $(OBJ3)
#
# link
# ----
#
trgl6: $(OBJB)
	f77 -o trgl6 $(OBJB)
trgl6_sqr: $(OBJA)
	f77 -o trgl6_sqr $(OBJA)
#
# all
# ---
#
all:
	make trgl6_sqr
	make trgl6
#
# compile
# -------
#
trgl6_sqr_dr.o: trgl6_sqr_dr.f
	f77 -c trgl6_sqr_dr.f
trgl6_disk.o: trgl6_disk.f
	f77 -c trgl6_disk.f
trgl6_sqr.o: trgl6_sqr.f
	f77 -c trgl6_sqr.f
trgl6.o: trgl6.f
	f77 -c trgl6.f
trgl6_sc.o: trgl6_sc.f
	f77 -c trgl6_sc.f
trgl6_ss.o: trgl6_ss.f
	f77 -c trgl6_ss.f
trgl6_sqr_geo.o: trgl6_sqr_geo.f
	f77 -c trgl6_sqr_geo.f
gauss_trgl.o: gauss_trgl.f
	f77 -c gauss_trgl.f
crvm_3d_2p.o: crvm_3d_2p.f
	f77 -c crvm_3d_2p.f
srf_int_3d_2p.o: srf_int_3d_2p.f
	f77 -c srf_int_3d_2p.f
#
# clean
# -----
#
clean:
	rm -f core
	rm -f $(OBJA) $(OBJB)
	rm -f trgl6_sqr trgl6 trgl6.stl trgl6.net
