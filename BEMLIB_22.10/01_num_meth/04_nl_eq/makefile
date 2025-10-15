#
# bairstow
# --------
#
bairstow: bairstow_dr.o bairstow.o bairstow_root.o
	f77 -o bairstow bairstow_dr.o bairstow.o bairstow_root.o
bairstow_dr.o: bairstow_dr.f
	f77 -c bairstow_dr.f
bairstow.o: bairstow.f
	f77 -c bairstow.f
bairstow_root.o: bairstow_root.f
	f77 -c bairstow_root.f
#
# cubic
# -----
#
cubic: cubic_dr.o cubic.o
	f77 -o cubic cubic_dr.o cubic.o
cubic_dr.o: cubic_dr.f
	f77 -c cubic_dr.f
cubic.o: cubic.f
	f77 -c cubic.f
#
#
# fpi1
# ----
#
fpi1: fpi1.f
	f77 -o fpi1 fpi1.f
#
#
# fpi2
# ----
#
fpi2: fpi2.f
	f77 -o fpi2 fpi2.f
#
#
# logistic
# --------
#
logistic: logistic.f
	f77 -o logistic logistic.f
#
#
# newton2 (two nonlinear equation)
# ------
#
OBJNEWTON2  = newton2_dr.o newton2.o newton2_fun.o 
newton2: $(OBJNEWTON2)
	f77 -o newton2 $(OBJNEWTON2)
newton2_dr.o: newton2_dr.f
	f77 -c newton2_dr.f
newton2.o: newton2.f
	f77 -c newton2.f
newton2_fun.o: newton2_fun.f
	f77 -c newton2_fun.f
#
#
# nol1 (one nonlinear equation)
# ----
#
OBJNOL1  = nol1.o newton1_2.o newton1_3.o muller.o secant.o
OBJNOL1A = newton1_fun.o muller_fun.o secant_fun.o
OBJNOL1B = bess_J0.o bess_J1.o
nol1: $(OBJNOL1) $(OBJNOL1A) $(OBJNOL1B)
	f77 -o nol1 $(OBJNOL1) $(OBJNOL1A) $(OBJNOL1B)
nol1.o: nol1.f
	f77 -c nol1.f
newton1_2.o: newton1_2.f
	f77 -c newton1_2.f
newton1_3.o: newton1_3.f
	f77 -c newton1_3.f
newton1_fun.o: newton1_fun.f
	f77 -c newton1_fun.f
muller.o: muller.f
	f77 -c muller.f
muller_fun.o: muller_fun.f
	f77 -c muller_fun.f
secant.o: secant.f
	f77 -c secant.f
secant_fun.o: secant_fun.f
	f77 -c secant_fun.f
bess_J0.o: bess_J0.f
	f77 -c bess_J0.f
bess_J1.o: bess_J1.f
	f77 -c bess_J1.f
#
#
# nol (arbitrary number of nonlinear equations)
# --
#
OBJNOL = nol.o nol_fun.o newton.o cramer_33.o newton_gel.o broyden.o 
nol: $(OBJNOL)
	f77 -o nol  $(OBJNOL)
broyden.o: broyden.f
	f77 -c broyden.f
cramer_33.o: cramer_33.f
	f77 -c cramer_33.f
newton.o: newton.f
	f77 -c newton.f
newton_gel: newton_gel.f
	f77 -c newton_gel.f 
nol.o: nol.f
	f77 -c nol.f
nol_fun.o: nol_fun.f
	f77 -c nol_fun.f
#
# quad
# ---------
#
quad: quad.f
	f77 -o quad quad.f
#
# quad_cc
# -------
#
quad_cc: quad.cc
	c++ -o quad_cc quad.cc
#
# quadc
# -----
#
quadc: quadc_dr.o quadc.o
	f77 -o quadc quadc_dr.o quadc.o
quadc_dr.o: quadc_dr.f
	f77 -c quadc_dr.f
quadc.o: quadc.f
	f77 -c quadc.f
#
#
# clean
# -----
#
clean:
	rm -f core nol1 nol1.out
	rm -f bairstow bairstow_dr.o bairstow.o bairstow_root.o 
	rm -f cubic_dr.o cubic.o cubic
	rm -f fpi1.o fpi1.out fpi1
	rm -f fpi2.o fpi2
	rm -f logistic.o logistic
	rm -f $(OBJNOL) nol nol.out
	rm -f $(OBJNOL1) $(OBJNOL1A) $(OBJNOL1B) nol1 nol1.out
	rm -f quad
	rm -f quadc quadc_dr.o quadc.o
	rm -f quad_cc
	rm -f newton2 $(OBJNEWTON2) newton2.out
#
#
# all fortran
# -----------
#
all:
	make bairstow
	make cubic
	make fpi1
	make fpi2
	make quad
	make quad_cc
	make quadc
	make logistic
	make nol
	make nol1
	make newton2
