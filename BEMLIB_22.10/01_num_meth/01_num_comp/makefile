#
# libraries
# ---------
#
LIBS = 
#
# compiler flags
# --------------
#
FLAGS =
#
#
# binary
# ------
#
binary: binary.f
	f77 -o binary binary.f $(LIBS)
#
# bits
# ----
#
bits: bits.f
	f77 -o bits bits.f
bits_cc: bits_cc.o powi.o
	g++ -o bits_cc bits_cc.o powi.o
powi.o: powi.cc
	g++ -c powi.cc
bits_cc.o: bits_cc.cc
	g++ -c bits_cc.cc
bits0_cc: bits0_cc.cc
	g++ -o bits0_cc bits0_cc.cc
#
# binomial_coefficient
# -------------
#
binomial_coeff: binomial_coeff.f
	f77 -o binomial_coeff binomial_coeff.f
binomial_coeff_cc: binomial_coeff.cc
	g++ -o binomial_coeff_cc binomial_coeff.cc
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
# cluster_2d
# ----------
#
OBJCL = cluster_2d_dr.o cluster_2d.o cluster_2d_link.o
cluster_2d: $(OBJCL)
	f77 -o cluster_2d $(OBJCL)
cluster_2d.o: cluster_2d.f
	f77 -c $(FLAGS) cluster_2d.f
cluster_2d_dr.o: cluster_2d_dr.f
	f77 -c $(FLAGS) cluster_2d_dr.f
cluster_2d_link.o: cluster_2d_link.f
	f77 -c $(FLAGS) cluster_2d_link.f
#
# gcd
# ---
#
gcd: gcd.f
	f77 -o gcd gcd.f
gcd_cc: gcd.cc
	g++ -o gcd_cc gcd.cc
#
# series
# ------
#
series: series.cc
	g++ -o series series.cc
#
# loss
# ---
#
loss: loss.cc
	g++ -o loss loss.cc
#
# golden ratio
# ------------
#
golden_ratio: golden_ratio.o 
	f77 -o golden_ratio golden_ratio.o 
#
# hello
# -----
#
hello_cc: hello.cc
	c++ -o hello_cc hello.cc
#
# mapping
# -------
#
mapping: mapping.f
	f77 -o mapping mapping.f
mapping_cc: mapping.cc
	c++ -o mapping_cc mapping.cc
#
#
# matrix_mult 
# -----------
#
matrix_mult: matrix_mult.f
	f77 -o matrix_mult matrix_mult.f
#
# pie
# ---
#
pie: pie.f
	f77 -o pie pie.f
#
#
# prime
# -----
#
prime: prime.f
	f77 -o prime prime.f
prime_cc: prime.cc
	c++ -o prime_cc prime.cc
#
# bsort
# -----
#
bsort: bsort.o random.o
	f77 -o bsort bsort.o random.o
random.o: random.f
	f77 -c random.f
#
#
# sum_fast
# --------
#
sum_fast: sum_fast.o 
	f77 -o sum_fast sum_fast.o 
#
# isort
# --------
#
isort: isort.o random.o
	f77 -o isort isort.o random.o
#
# richardson
# ----------
#
richardson: richardson.o 
	f77 -o richardson richardson.o 
#
#
# von_koch
# ----------
#
OBJVK = von_koch_dr.o von_koch.o
von_koch: $(OBJVK)
	f77 -o von_koch $(OBJVK)
von_koch_dr.o: von_koch_dr.f
	f77 -c von_koch_dr.f
von_koch.o: von_koch.f
	f77 -c von_koch.f
#
# sum the elements of a matrix
# ----------------------------
#
msum: msum.f
	f77 -o msum msum.f
#
# all fortran
# -----------
#
all:
	make bits
	make bits_cc
	make bits0_cc
	make binary 
	make cluster_2d
	make binomial_coeff
	make binomial_coeff_cc
	make gcd
	make golden_ratio
	make matrix_mult
	make mapping
	make pie
	make prime
	make bsort
	make isort
	make richardson
	make sum_fast
	make von_koch
	make msum
#
#
# all c++
# -------
#
#
all_cc:
	make bits_cc
	make bits0_cc
	make combinatorial_cc
	make mapping_cc
	make hello_cc
	make prime_cc
	make series
#
# clean
# -----
#
clean:
	rm -f bits.o bits bits_cc.o bits_cc powi.o bits0_cc.o bits0_cc
	rm -f binary.o binary 
	rm -f $(OBJCL) cluster_2d
	rm -f binomial_coeff
	rm -f binomial_coeff_cc 
	rm -f gcd.o gcd gcd_cc 
	rm -f sum_fast.o sum_fast 
	rm -f golden_ratio.o golden_ratio
	rm -f hello_cc
	rm -f matrix_mult.o matrix_mult
	rm -f mapping.o mapping mapping_cc
	rm -f pie.o pie
	rm -f prime.o prime prime_cc
	rm -f bsort.o bsort random.o
	rm -f isort.o isort random.o
	rm -f richardson.o richardson
	rm -f $(OBJVK) von_koch
	rm -f msum
	rm -f series
	rm -f loss
	rm -f hello
