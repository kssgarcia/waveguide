OBJ1 = prtcl_ax.o precondition.o reduce.o
OBJ2 = prtcl_ax_geo.o prtcl_ax_slp.o 
OBJ3 = sgf_ax_fs.o sgf_ax_w.o sgf_ax_1p_ct.o
OBJ4 = ell_int.o gauss_leg.o bess_I01K01.o
OBJ5 = gel_inv.o
OBJ6 = prtcl_ax_sys.o prtcl_ax_str.o prtcl_ax_vel.o 
OBJ  = $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6)
#
# link
# ----
#
prtcl_ax: $(OBJ)
	f77 -o prtcl_ax $(OBJ)
#
# compile
# -------
#
prtcl_ax.o:prtcl_ax.f
	f77 -c prtcl_ax.f
prtcl_ax_geo.o: prtcl_ax_geo.f
	f77 -c prtcl_ax_geo.f
prtcl_ax_sys.o:prtcl_ax_sys.f
	f77 -c prtcl_ax_sys.f
precondition.o: precondition.f
	f77 -c precondition.f
reduce.o: reduce.f
	f77 -c reduce.f
prtcl_ax_slp.o: prtcl_ax_slp.f
	f77 -c prtcl_ax_slp.f
prtcl_ax_str.o: prtcl_ax_str.f
	f77 -c prtcl_ax_str.f
prtcl_ax_vel.o: prtcl_ax_vel.f
	f77 -c prtcl_ax_vel.f
gauss_leg.o: gauss_leg.f
	f77 -c gauss_leg.f
gel_inv.o: gel_inv.f
	f77 -c gel_inv.f
sgf_ax_fs.o: sgf_ax_fs.f
	f77 -c sgf_ax_fs.f
sgf_ax_w.o: sgf_ax_w.f
	f77 -c sgf_ax_w.f
sgf_ax_1p_ct.o: sgf_ax_1p_ct.f
	f77 -c sgf_ax_1p_ct.f
ell_int.o: ell_int.f
	f77 -c ell_int.f
bess_I01K01.o: bess_I01K01.f
	f77 -c bess_I01K01.f
#
# clean
# ----
#
clean:
	rm -f core PLOTDAT
	rm -f $(OBJ) prtcl_ax
	rm -f PLOTDAT
	rm -f prtcl_ax.str prtcl_ax.trc prtcl_ax.error
	rm -f prtcl_ax.xy
#
# purge
# ----
#
purge:
	rm -f core PLOTDAT AAA-xmgr
	rm -f $(OBJ)
	rm -f PLOTDAT
	rm -f prtcl_ax.str prtcl_ax.trc prtcl_ax.error
	rm -f prtcl_ax.xy
#
# all
# ----
#
all:
	make prtcl_ax
#
# tar
# ---
#
tar:
	tar cvf prtcl_ax.tar makefile *.f *.dat

