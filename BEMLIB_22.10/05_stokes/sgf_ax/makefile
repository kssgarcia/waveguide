#
# sgf_ax_fs
# ---------
#
OBJ1 = sgf_ax_fs_dr.o sgf_ax_fs.o ell_int.o
sgf_ax_fs: $(OBJ1)
	f77 -o sgf_ax_fs $(OBJ1)
sgf_ax_fs_dr.o: sgf_ax_fs_dr.f
	f77 -c sgf_ax_fs_dr.f
sgf_ax_fs.o: sgf_ax_fs.f
	f77 -c sgf_ax_fs.f
#
# sgf_ax_w
# --------
#
OBJ2 = sgf_ax_w_dr.o sgf_ax_w.o ell_int.o sgf_ax_fs.o
sgf_ax_w: $(OBJ2)
	f77 -o sgf_ax_w $(OBJ2)
sgf_ax_w_dr.o: sgf_ax_w_dr.f
	f77 -c sgf_ax_w_dr.f
sgf_ax_w.o: sgf_ax_w.f
	f77 -c sgf_ax_w.f
#
# sgf_ax_1p
# ---------
#
OBJ3 = sgf_ax_1p_dr.o sgf_ax_1p.o sgf_ax_fs.o ell_int.o 
sgf_ax_1p: $(OBJ3)
	f77 -o sgf_ax_1p $(OBJ3)
sgf_ax_1p_dr.o: sgf_ax_1p_dr.f
	f77 -c sgf_ax_1p_dr.f
sgf_ax_1p.o: sgf_ax_1p.f
	f77 -c sgf_ax_1p.f
#
# sgf_ax_1p_ct
# ------------
#
OBJ4 = sgf_ax_1p_ct_dr.o sgf_ax_1p_ct.o ell_int.o bess_I01K01.o sgf_ax_fs.o
sgf_ax_1p_ct: $(OBJ4)
	f77 -o sgf_ax_1p_ct $(OBJ4)
sgf_ax_1p_ct_dr.o: sgf_ax_1p_ct_dr.f
	f77 -c sgf_ax_1p_ct_dr.f
sgf_ax_1p_ct.o: sgf_ax_1p_ct.f
	f77 -c sgf_ax_1p_ct.f
bess_I01K01.o:
	f77 -c bess_I01K01.f
#
# ell_int
# -------
#
ell_int.o:
	f77 -c ell_int.f
#
# all
# ---
#
all:
	make sgf_ax_fs
	make sgf_ax_w
	make sgf_ax_1p
	make sgf_ax_1p_ct
#
# clean
# -----
#
clean:
	rm -f core PLOTDAT
	rm -f $(OBJ1) sgf_ax_fs
	rm -f $(OBJ2) sgf_ax_w
	rm -f $(OBJ3) sgf_ax_1p
	rm -f $(OBJ4) sgf_ax_1p_ct
#
# purge
# -----
#
purge:
	rm -f core PLOTDAT
	rm -f $(OBJ1)
	rm -f $(OBJ2)
	rm -f $(OBJ3)
	rm -f $(OBJ4)
#
# tar
# -----
#
tar:
	tar cvf sgf_ax.tar makefile *.f
