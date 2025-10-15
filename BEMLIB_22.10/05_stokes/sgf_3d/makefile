#
# sgf_3d_fs
# ---------
#
OBJ1= sgf_3d_fs_dr.o sgf_3d_fs.o gauss_sph.o
sgf_3d_fs: $(OBJ1)
	f77 -o sgf_3d_fs $(OBJ1)
sgf_3d_fs_dr.o: sgf_3d_fs_dr.f
	f77 -c sgf_3d_fs_dr.f
sgf_3d_fs.o: sgf_3d_fs.f
	f77 -c sgf_3d_fs.f
#
# sgf_3d_w
# ---------
#
OBJ2 = sgf_3d_w_dr.o sgf_3d_w.o gauss_sph.o
sgf_3d_w: $(OBJ2)
	f77 -o sgf_3d_w $(OBJ2)
sgf_3d_w_dr.o: sgf_3d_w_dr.f
	f77 -c sgf_3d_w_dr.f
sgf_3d_w.o: sgf_3d_w.f
	f77 -c sgf_3d_w.f
#
# sgf_3d_2p_w
# -----------
#
OBJ_3d_w = sgf_3d_2p_w_dr.o sgf_3d_2p_w.o sgf_3d_w.o gauss_sph.o

sgf_3d_2p_w: $(OBJ_3d_w)
	f77 -o sgf_3d_2p_w $(OBJ_3d_w)
sgf_3d_2p_w_dr.o: sgf_3d_2p_w_dr.f
	f77 -c sgf_3d_2p_w_dr.f
sgf_3d_2p_w.o: sgf_3d_2p_w.f
	f77 -c sgf_3d_2p_w.f
#
# sgf_3d_fs_sph
# -------------
#
OBJ3= sgf_3d_sph_dr.o sgf_3d_sph.o gauss_sph.o
sgf_3d_sph: $(OBJ3) 
	f77 -o sgf_3d_sph $(OBJ3)
sgf_3d_sph_dr.o: sgf_3d_sph_dr.f
	f77 -c sgf_3d_sph_dr.f
sgf_3d_sph.o: sgf_3d_sph.f
	f77 -c sgf_3d_sph.f
#
# sgf_3d_2p
# ---------
#
OBJ_3d_2p = sgf_3d_2p_dr.o sgf_3d_2p.o sgf_3d_2p_aux.o sgf_3d_2p_ewald.o

sgf_3d_2p: $(OBJ_3d_2p) 
	f77 -o sgf_3d_2p $(OBJ_3d_2p)
sgf_3d_2p_aux: sgf_3d_2p_aux.f
	f77 -c sgf_3d_2p_aux.f
sgf_3d_2p_dr.o: sgf_3d_2p_dr.f
	f77 -c sgf_3d_2p_dr.f
sgf_3d_2p.o: sgf_3d_2p.f
	f77 -c sgf_3d_2p.f
sgf_3d_2p_ewald.o: sgf_3d_2p_ewald.f
	f77 -c sgf_3d_2p_ewald.f
#
# sgf_3d_3p
# ---------
#
OBJ_3d_3p = sgf_3d_3p_dr.o sgf_3d_3p.o sgf_3d_3p_ewald.o sgf_3d_3p_qqq.o sgf_3d_3p_vvv.o

sgf_3d_3p: $(OBJ_3d_3p) gauss_sph.o
	f77 -o sgf_3d_3p $(OBJ_3d_3p) gauss_sph.o
sgf_3d_3p_dr.o: sgf_3d_3p_dr.f
	f77 -c sgf_3d_3p_dr.f
sgf_3d_3p.o: sgf_3d_3p.f
	f77 -c sgf_3d_3p.f
sgf_3d_3p_ewald.o: sgf_3d_3p_ewald.f
	f77 -c sgf_3d_3p_ewald.f
sgf_3d_3p_qqq.o: sgf_3d_3p_qqq.f
	f77 -c sgf_3d_3p_qqq.f
sgf_3d_3p_vvv.o: sgf_3d_3p_vvv.f
	f77 -c sgf_3d_3p_vvv.f
#
# Gauss quadrature
# ----------------
#
gauss_sph.o:
	f77 -c gauss_sph.f
#
# all
# ---
#
all:
	make sgf_3d_fs
	make sgf_3d_w
	make sgf_3d_sph
	make sgf_3d_2p
	make sgf_3d_3p
	make sgf_3d_2p_w
#
# tar
# ---
#
tar:
	tar cvf sgf_3d.tar makefile *.f *.dat
#
# clean
# -----
#
clean:
	rm -f core
	rm -f $(OBJ1) sgf_3d_fs
	rm -f $(OBJ2) sgf_3d_w
	rm -f $(OBJ3) sgf_3d_sph
	rm -f $(OBJ_3d_2p) sgf_3d_2p
	rm -f $(OBJ_3d_3p) sgf_3d_3p
	rm -f $(OBJ_3d_w) sgf_3d_2p_w
	rm -f sgf.out
#
# purge
# -----
#
purge:
	rm -f core
	rm -f $(OBJ1) 
	rm -f $(OBJ2) 
	rm -f $(OBJ3) 
	rm -f $(OBJ_3d_2p) 
	rm -f $(OBJ_3d_3p)
	rm -f $(OBJ_3d_w)
	rm -f sgf.out
