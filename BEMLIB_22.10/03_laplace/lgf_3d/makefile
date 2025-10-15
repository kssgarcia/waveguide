#
# lgf_3d_fs
# ---------
#
OBJ1= lgf_3d_fs_dr.o lgf_3d_fs.o gauss_sph.o
lgf_3d_fs: $(OBJ1) 
	f77 -o lgf_3d_fs $(OBJ1)
lgf_3d_fs_dr.o: lgf_3d_fs_dr.f
	f77 -c lgf_3d_fs_dr.f
lgf_3d_fs.o: lgf_3d_fs.f
	f77 -c lgf_3d_fs.f
gauss_sph.o: gauss_sph.f
	f77 -c gauss_sph.f
#
#
# lgf_3d_w
# --------
#
OBJ2= lgf_3d_w_dr.o lgf_3d_w.o gauss_sph.o
lgf_3d_w: $(OBJ2)
	f77 -o lgf_3d_w $(OBJ2)
lgf_3d_w_dr.o: lgf_3d_w_dr.f
	f77 -c lgf_3d_w_dr.f
lgf_3d_w.o: lgf_3d_w.f
	f77 -c lgf_3d_w.f
#
# lgf_3d_sph
# ----------
#
OBJ4= lgf_3d_sph_dr.o lgf_3d_sph.o gauss_sph.o
lgf_3d_sph: $(OBJ4)
	f77 -o lgf_3d_sph $(OBJ4)
lgf_3d_sph_dr.o: lgf_3d_sph_dr.f
	f77 -c lgf_3d_sph_dr.f
lgf_3d_sphere.o: lgf_3d_sphere.f
	f77 -c lgf_3d_sphere.f
#
# lgf_3d_2p
# ----------
#
OBJ5 = lgf_3d_2p_dr.o lgf_3d_2p.o lgf_3d_2p_ewald.o gauss_sph.o 
lgf_3d_2p: $(OBJ5)
	f77 -o lgf_3d_2p $(OBJ5)
lgf_3d_2p_dr.o: lgf_3d_2p_dr.f
	f77 -c lgf_3d_2p_dr.f
lgf_3d_2p.o: lgf_3d_2p.f
	f77 -c lgf_3d_2p.f
lgf_3d_2p_ewald.o: lgf_3d_2p_ewald.f
	f77 -c lgf_3d_2p_ewald.f
#
#
#
# lgf_3d_2p_glut
# --------------
#
OBJ6= lgf_3d_2p_glut.o lgf_3d_2p.o lgf_3d_2p_int.o lgf_3d_2p_ewald.o
lgf_3d_2p_glut: $(OBJ6)
	f77 -o lgf_3d_2p_glut $(OBJ6)
lgf_3d_2p_glut.o: lgf_3d_2p_glut.f
	f77 -c lgf_3d_2p_glut.f
lgf_3d_2p_int.o: lgf_3d_2p_int.f
	f77 -c lgf_3d_2p_int.f
#
# clean
# -----
#
clean:
	rm -f core PLOTDAT 
	rm -f $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6)
	rm -f lgf_3d_fs lgf_3d_w lgf_3d_sph lgf_3d_2p
	rm -f lgf_3d_2p_glut
#
# purge
# -----
#
purge:
	rm -f core PLOTDAT 
	rm -f $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6)
#
# all
# ---
#
all:
	make lgf_3d_fs
	make lgf_3d_w
	make lgf_3d_sph
	make lgf_3d_2p
	make lgf_3d_2p_glut
#
# tar
# ---
#
tar:
	tar cvf lgf_3d.tar makefile *.f
