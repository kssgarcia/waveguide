all:	
	make chimney_s0
	make chimney_s1
	make chimney_u1
#
# chimneys
# --------
#
chimney_s0: chimney_s0.o
	f77 -o chimney_s0 chimney_s0.o
chimney_s1: chimney_s1.o
	f77 -o chimney_s1 chimney_s1.o
chimney_u1: chimney_u1.o
	f77 -o chimney_u1 chimney_u1.o
#
# clean
# ---
#
clean:
	rm -f chimney_s0.o chimney_s0.out chimney_s0
	rm -f chimney_s1.o chimney_s1.out chimney_s1
	rm -f chimney_u1.o chimney_u1.out1 chimney_u1.out2 chimney_u1
