#
# lgf_ax_fs
# ---------
#
OBJfs = lgf_ax_fs_dr.o lgf_ax_fs.o ell_int.o
lgf_ax_fs: $(OBJfs) 
	f77 -o lgf_ax_fs $(OBJfs)
lgf_ax_fs_dr.o: lgf_ax_fs_dr.f
	f77 -c lgf_ax_fs_dr.f
lgf_ax_fs.o: lgf_ax_fs.f
	f77 -c lgf_ax_fs.f
#
# lgf_ax_w
# ---------
#
OBJw = lgf_ax_w_dr.o lgf_ax_w.o ell_int.o
lgf_ax_w: $(OBJw) 
	f77 -o lgf_ax_w $(OBJw)
lgf_ax_w_dr.o:lgf_ax_w_dr.f
	f77 -c lgf_ax_w_dr.f
lgf_ax_w.o: lgf_ax_w.f
	f77 -c lgf_ax_w.f
#
# common
# ------
# 
ell_int.o: ell_int.f
	f77 -c ell_int.f
#
# clean
# -----
# 
clean:
	rm -f core
	rm -f $(OBJfs) $(OBJw)
	rm -f lgf_ax_fs
	rm -f lgf_ax_w
#
# purge:
# -----
# 
purge:
	rm -f core
	rm -f $(OBJfs) $(OBJw)
#
# all:
# ---
#
all:
	make lgf_ax_fs
	make lgf_ax_w
#
# tar
# ---
#
tar:
	tar cvf lgf_ax.tar makefile *.f
