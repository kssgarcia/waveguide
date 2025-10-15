#
# arc_2d
# ------
#
O_arc_2d = arc_2d_dr.o arc_2d.o
arc_2d: $(O_arc_2d)
	f77 -o arc_2d $(O_arc_2d)
arc_2d_dr.o: arc_2d_dr.f
	f77 -c arc_2d_dr.f
arc_2d.o: arc_2d.f
	f77 -c arc_2d.f
#
# arc_3d
# ------
#
O_arc_3d = arc_3d_dr.o arc_3d.o cramer_33.o
arc_3d: $(O_arc_3d)
	f77 -o arc_3d $(O_arc_3d)
arc_3d_dr.o: arc_3d_dr.f
	f77 -c arc_3d_dr.f
arc_3d.o: arc_3d.f
	f77 -c arc_3d.f
cramer_33.o: cramer_33.f
	f77 -c cramer_33.f
#
# arc_ax
# ------
#
#
O_arc_ax = arc_ax_dr.o arc_ax.o
arc_ax: $(O_arc_ax)
	f77 -o arc_ax $(O_arc_ax)
arc_ax_dr.o: arc_ax_dr.f
	f77 -c arc_ax_dr.f
arc_ax.o: arc_ax.f
	f77 -c arc_ax.f
#
#
# bih_2d
# ------
#
bih_2d: bih_2d.o lapl_2d.o
	f77 -o bih_2d bih_2d.o lapl_2d.o 
bih_2d.o: bih_2d.f 
	f77 -c bih_2d.f 
#
# cheby
# -----
#
cheb: cheb.o
	f77 -o cheb cheb.o 
#
# hermite
# -------
#
hermite: hermite.o hermite_card.o hermite_table.o
	f77 -o hermite hermite.o hermite_card.o hermite_table.o
hermite_card.o: hermite_card.f 
	f77 -c hermite_card.f 
hermite_table.o: hermite_table.f 
	f77 -c hermite_table.f 
#
# interpolate
# -----------
#
interpolate_o = interpolate.o lagrange.o neville.o aitken.o newton.o fornberg.o
interpolate: $(interpolate_o) 
	f77 -o interpolate $(interpolate_o)
lagrange.o: lagrange.f 
	f77 -c lagrange.f 
neville.o: neville.f 
	f77 -c neville.f 
aitken.o: aitken.f 
	f77 -c aitken.f 
newton.o: newton.f 
	f77 -c newton.f 
fornberg.o: fornberg.f 
	f77 -c fornberg.f 
#
# lagrange_int
#  -----------
#
lagrange_int: lagrange_int.o lobatto.o
	f77 -o lagrange_int lagrange_int.f lobatto.o
lobatto.o: lobatto.f 
	f77 -c lobatto.f 
#
# lagrange_es
# -----------
#
lagrange_es: lagrange_es.f
	f77 -o lagrange_es lagrange_es.f
#
#
# lapl_2d
# -------
#
lapl_2d: lapl_2d_dr.o lapl_2d.o
	f77 -o lapl_2d lapl_2d_dr.o lapl_2d.o 
lapl_2d_dr.o: lapl_2d_dr.f 
	f77 -c lapl_2d_dr.f 
lapl_2d.o: lapl_2d.f 
	f77 -c lapl_2d.f 
#
# cubic splines
# -------------
#
OBJ_splc1 = splc_dr.o splc_pr.o splc_nt.o splc_clm.o 
OBJ_splc2 = thomas.o thomas_pr.o 
OBJ_splc9 = splc_pr_dr.o splc_pr.o thomas_pr.o 
OBJ_splc  = $(OBJ_splc1) $(OBJ_splc2)

splc: $(OBJ_splc)
	f77 -o splc $(OBJ_splc)
splc_pr_dr: $(OBJ_splc9)
	f77 -o splc_pr_dr $(OBJ_splc9)
splc_pr_dr.o: splc_pr_dr.f 
	f77 -c splc_pr_dr.f 

splc_dr.o: splc_dr.f 
	f77 -c splc_dr.f 
splc_pr.o: splc_pr.f 
	f77 -c splc_pr.f 
splc_nt.o: splc_nt.f 
	f77 -c splc_nt.f 
splc_clm.o: splc_clm.f 
	f77 -c splc_clm.f 
thomas.o: thomas.f 
	f77 -c thomas.f 
thomas_pr.o: thomas_pr.f 
	f77 -c thomas_pr.f 
#
# newton_cc
# ---------
#
newton_cc: newton_cc.cc
	g++ -o newton_cc newton_cc.cc
#
# newton_c
# ---------
#
newton_c: newton_c.c
	g++ -o newton_c newton_c.c
#
# all
# ---
#
all:
	make arc_2d
	make arc_3d
	make arc_ax
	make bih_2d
	make cheb
	make lagrange_int
	make lagrange_es
	make lapl_2d
	make interpolate
	make hermite
	make splc
	make newton_cc
	make newton_c
#
# clean
# -----
#
clean:
	rm -f core PLOTDAT AAA-xmgr
	rm -f $(O_arc_2d) arc_2d
	rm -f $(O_arc_3d) arc_3d
	rm -f $(O_arc_ax) arc_ax
	rm -f bih_2d.o lapl_2d.o bih_2d
	rm -f cheb.o cheb
	rm -f hermite.o hermite_card.o hermite_table.o hermite
	rm -f $(interpolate_o) interpolate
	rm -f lagrange_int.o lobatto.o lagrange_int
	rm -f lagrange_es.o lagrange_es
	rm -f lapl_2d.o lapl_2d_dr.o lapl_2d
	rm -f $(OBJ_splc) splc.out splc
	rm -f $(OBJ_splc9) splc_periodic
	rm -f newton_cc newton_c
#
# tar
# ---
#
tar:
	tar cvf 06_interp_diff.tar makefile *.f *dat 06*
