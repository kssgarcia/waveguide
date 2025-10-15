#
# lgf_2d_fs
# ----------
#
#
OBJfs = lgf_2d_fs_dr.o lgf_2d_fs.o 
lgf_2d_fs:  $(OBJfs) 
	f77 -o lgf_2d_fs $(OBJfs)
lgf_2d_fs_dr.o: lgf_2d_fs_dr.f
	f77 -c lgf_2d_fs_dr.f
lgf_2d_fs.o: lgf_2d_fs.f
	f77 -c lgf_2d_fs.f
#
# lgf_2d_1p
# --------
#
OBJ1p = lgf_2d_1p_dr.o lgf_2d_1p.o 
#
lgf_2d_1p: $(OBJ1p) 
	f77 -o lgf_2d_1p $(OBJ1p)
lgf_2d_1p_dr.o: lgf_2d_1p_dr.f
	f77 -c lgf_2d_1p_dr.f
lgf_2d_1p.o: lgf_2d_1p.f
	f77 -c lgf_2d_1p.f
#
#
# lgf_2d_w
# --------
#
OBJw = lgf_2d_w_dr.o lgf_2d_w.o 
#
lgf_2d_w: $(OBJw) 
	f77 -o lgf_2d_w $(OBJw)
lgf_2d_w_dr.o: lgf_2d_w_dr.f
	f77 -c lgf_2d_w_dr.f
lgf_2d_w.o: lgf_2d_w.f
	f77 -c lgf_2d_w.f
#
# lgf_2d_1p_w
# -----------
#
OBJ1pw = lgf_2d_1p_w_dr.o lgf_2d_1p_w.o 
#
lgf_2d_1p_w: $(OBJ1pw) 
	f77 -o lgf_2d_1p_w $(OBJ1pw)
lgf_2d_1p_w_dr.o:
	f77 -c lgf_2d_1p_w_dr.f
lgf_2d_1p_w.o:
	f77 -c lgf_2d_1p_w.f
#
#
# lgf_2d_ww
# ---------
#
OBJww = lgf_2d_ww_dr.o lgf_2d_ww.o 
#
lgf_2d_ww: $(OBJww) 
	f77 -o lgf_2d_ww $(OBJww)
lgf_2d_ww_dr.o: lgf_2d_ww_dr.f
	f77 -c lgf_2d_ww_dr.f
lgf_2d_ww.o: lgf_2d_ww.f
	f77 -c lgf_2d_ww.f
#
#
# lgf_2d_www
# ----------
#
OBJwww = lgf_2d_www_dr.o lgf_2d_www.o 
#
lgf_2d_www: $(OBJwww) 
	f77 -o lgf_2d_www $(OBJwww)
lgf_2d_www_dr.o:
	f77 -c lgf_2d_www_dr.f
lgf_2d_www.o:
	f77 -c lgf_2d_www.f
#
#
# lgf_2d_crc
# ----------
#
OBJcrc = lgf_2d_crc_dr.o lgf_2d_crc.o 
#
lgf_2d_crc: $(OBJcrc) 
	f77 -o lgf_2d_crc $(OBJcrc)
lgf_2d_crc_dr.o:lgf_2d_crc_dr.f
	f77 -c lgf_2d_crc_dr.f
lgf_2d_crc.o: lgf_2d_crc.f
	f77 -c lgf_2d_crc.f
#
# clean
# -----
#
clean:
	rm -f core
	rm -f $(OBJfs)  lgf_2d_fs
	rm -f $(OBJ1p)  lgf_2d_1p
	rm -f $(OBJw)   lgf_2d_w
	rm -f $(OBJ1pw) lgf_2d_1p_w
	rm -f $(OBJww)  lgf_2d_ww
	rm -f $(OBJwww) lgf_2d_www
	rm -f $(OBJcrc) lgf_2d_crc
#
# purge
# -----
#
purge:
	rm -f core
	rm -f $(OBJfs)
	rm -f $(OBJ1p)
	rm -f $(OBJw)
	rm -f $(OBJ1pw)
	rm -f $(OBJww)
	rm -f $(OBJwww)
	rm -f $(OBJcrc)
#
# all:
# ---
#
all:
	make lgf_2d_fs
	make lgf_2d_1p
	make lgf_2d_w
	make lgf_2d_1p_w
	make lgf_2d_ww
	make lgf_2d_www
	make lgf_2d_crc
#
# tar
# ---
#
tar:
	tar cvf lgf_2d.tar makefile *.f
