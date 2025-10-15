#
# Green's functions of 2D Stokes flow
#
#
# sgf_2d_fs
# ---------
#
OBJ1 = sgf_2d_fs_dr.o sgf_2d_fs.o 
sgf_2d_fs: $(OBJ1) 
	f77 -o sgf_2d_fs $(OBJ1)
sgf_2d_fs_dr.o: sgf_2d_fs_dr.f
	f77 -c sgf_2d_fs_dr.f
sgf_2d_fs.o: sgf_2d_fs.f
	f77 -c sgf_2d_fs.f
#
# sgf_2d_w
# --------
#
OBJ2= sgf_2d_w_dr.o sgf_2d_w.o 
sgf_2d_w: $(OBJ2)
	f77 -o sgf_2d_w $(OBJ2)
sgf_2d_w_dr.o: sgf_2d_w_dr.f
	f77 -c sgf_2d_w_dr.f
sgf_2d_w.o: sgf_2d_w.f
	f77 -c sgf_2d_w.f
#
# sgf_2d_1p
# ---------
#
OBJ3 = sgf_2d_1p_dr.o sgf_2d_1p.o 
sgf_2d_1p: $(OBJ3) 
	f77 -o sgf_2d_1p $(OBJ3)
sgf_2d_1p_dr.o: sgf_2d_1p_dr.f
	f77 -c sgf_2d_1p_dr.f
sgf_2d_1p.o: sgf_2d_1p.f
	f77 -c sgf_2d_1p.f
#
# sgf_2d_1p_w
# -----------
#
OBJ4 = sgf_2d_1p_w_dr.o sgf_2d_1p_w.o 
sgf_2d_1p_w: $(OBJ4) 
	f77 -o sgf_2d_1p_w $(OBJ4)
sgf_2d_1p_w_dr.o: sgf_2d_1p_w_dr.f
	f77 -c sgf_2d_1p_w_dr.f
sgf_2d_1p_w.o: sgf_2d_1p_w.f
	f77 -c sgf_2d_1p_w.f
#
# sgf_2d_1p_ww
# ------------
#
OBJ5A = sgf_2d_1p_ww_dr.o sgf_2d_1p_ww.o
OBJ5B = sgf_2d_1p_ww_fga.o sgf_2d_1p_ww_int.o sgf_2d_1p_w.o 
OBJ5  = $(OBJ5A) $(OBJ5B) 
sgf_2d_1p_ww: $(OBJ5)
	f77 -o sgf_2d_1p_ww $(OBJ5)
sgf_2d_1p_ww_dr.o: sgf_2d_1p_ww_dr.f
	f77 -c sgf_2d_1p_ww_dr.f
sgf_2d_1p_ww.o: sgf_2d_1p_ww.f
	f77 -c sgf_2d_1p_ww.f
sgf_2d_1p_ww_fga.o: sgf_2d_1p_ww_fga.f
	f77 -c sgf_2d_1p_ww_fga.f
sgf_2d_1p_ww_int.o: sgf_2d_1p_ww_int.f
	f77 -c sgf_2d_1p_ww_int.f
#
# sgf_2d_1p_ww_glut
# -----------------
#
OBJ6A = sgf_2d_1p_ww_glut.o sgf_2d_1p_ww.o 
OBJ6B = sgf_2d_1p_ww_fga.o sgf_2d_1p_ww_int.o sgf_2d_1p_w.o 
OBJ6  = $(OBJ6A) $(OBJ6B) 
sgf_2d_1p_ww_glut: $(OBJ6)
	f77 -o sgf_2d_1p_ww_glut $(OBJ6)
sgf_2d_1p_ww_glut.o: sgf_2d_1p_ww_glut.f
	f77 -c sgf_2d_1p_ww_glut.f
#
# sgf_2d_2p
# ---------
#
OBJ7A = sgf_2d_2p_dr.o sgf_2d_2p.o sgf_2d_2p_int.o exp_int.o
OBJ7B = sgf_2d_2p_ewald.o sgf_2d_2p_qqq.o sgf_2d_2p_vvv.o
OBJ7  = $(OBJ7A) $(OBJ7B) 
sgf_2d_2p: $(OBJ7)
	f77 -o sgf_2d_2p $(OBJ7)
sgf_2d_2p_dr.o: sgf_2d_2p_dr.f
	f77 -c sgf_2d_2p_dr.f
sgf_2d_2p.o: sgf_2d_2p.f
	f77 -c sgf_2d_2p.f
sgf_2d_2p_ewald.o: sgf_2d_2p_ewald.f
	f77 -c sgf_2d_2p_ewald.f
sgf_2d_2p_vvv.o: sgf_2d_2p_vvv.f
	f77 -c sgf_2d_2p_vvv.f
sgf_2d_2p_qqq.o: sgf_2d_2p_qqq.f
	f77 -c sgf_2d_2p_qqq.f
exp_int.o: exp_int.f
	f77 -c exp_int.f
sgf_2d_2p_int.o: sgf_2d_2p_int.f
	f77 -c sgf_2d_2p_int.f
#
# sgf_2d_2p_glut
# --------------
#
OBJ8A = sgf_2d_2p_glut.o sgf_2d_2p.o sgf_2d_2p_int.o exp_int.o
OBJ8B = sgf_2d_2p_ewald.o sgf_2d_2p_qqq.o sgf_2d_2p_vvv.o
OBJ8  = $(OBJ8A) $(OBJ8B) 
sgf_2d_2p_glut: $(OBJ8)
	f77 -o sgf_2d_2p_glut $(OBJ8)
sgf_2d_2p_glut.o: sgf_2d_2p_glut.f
	f77 -c sgf_2d_2p_glut.f
#
# all
# ---
#
all:
	make sgf_2d_fs
	make sgf_2d_w
	make sgf_2d_1p
	make sgf_2d_1p_w
	make sgf_2d_1p_ww
	make sgf_2d_2p
	make sgf_2d_2p_glut
	make sgf_2d_1p_ww_glut
#
# clean
# -----
#
clean:
	rm -f core
	rm -f $(OBJ1) sgf_2d_fs
	rm -f $(OBJ2) sgf_2d_w
	rm -f $(OBJ3) sgf_2d_1p
	rm -f $(OBJ4) sgf_2d_1p_w
	rm -f $(OBJ5) sgf_2d_1p_ww sgf_2d_1p_ww.out
	rm -f $(OBJ6) sgf_2d_1p_ww_glut
	rm -f $(OBJ7) sgf_2d_2p
	rm -f $(OBJ8) sgf_2d_2p_glut
#
# purge
# -----
#
purge:
	rm -f core
	rm -f $(OBJ1)
	rm -f $(OBJ2)
	rm -f $(OBJ3)
	rm -f $(OBJ4)
	rm -f $(OBJ5) sgf_2d_1p_ww.out
	rm -f $(OBJ6)
	rm -f $(OBJ7)
	rm -f $(OBJ8)
#
# tar
# ---
#
tar:
	tar cvf sgf_2d.tar makefile *.f *.dat
