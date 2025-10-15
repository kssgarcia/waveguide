#
# object files
#
OBJ1 = em_2d.o em_2d_vel.o em_2d_dlp.o taylor.o
OBJ2 = prd_2d.o prd_2d_interp.o arc_2d.o gel_large.o gel_small.o
OBJ3 = sgf_2d_fs.o sgf_2d_1p_w.o sgf_2d_1p_ww.o sgf_2d_1p_ww_fga.o 
OBJ4 = sgf_2d_2p_int.o sgf_2d_w.o
OBJ5 = gauss_leg.o exp_int.o verbal.o
OBJ6 = em_2d_rhe.o em_2d_pdfl.o deflation.o
OBJ7 = sgf_2d_2p.o sgf_2d_2p_ewald.o sgf_2d_2p_qqq.o sgf_2d_2p_vvv.o
OBJ8 = em_2d_slp.o slp_arc_host.o slp_arc_blended.o
OBJ9 = em_2d_cd.o interp_en.o interp_mn.o
OBJ  = $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) $(OBJ7) $(OBJ8) $(OBJ9)
#
#
# libraries
# ---------
#
LIBS = -lm 
#
#
# compiler flags
# --------------
#
FLAGS =
#
# link
# ----
#
em_2d: $(OBJ)
	f77 -o em_2d $(OBJ) $(LIBS)
#
# compile
# -------
#
em_2d.o: em_2d.f
	f77 -c $(flags) em_2d.f
verbalo: verbal.f
	f77 -c $(flags) verbal.f
arc_2d.o: arc_2d.f
	f77 -c $(flags) arc_2d.f
em_2d_cd.o: em_2d_cd.f
	f77 -c $(flags) em_2d_cd.f
interp_en.o: interp_en.f
	f77 -c $(flags) interp_en.f
interp_mn.o: interp_mn.f
	f77 -c $(flags) interp_mn.f
em_2d_vel.o: em_2d_vel.f
	f77 -c $(flags) em_2d_vel.f
deflation.o: deflation.f
	f77 -c $(flags) deflation.f
em_2d_slp.o: em_2d_slp.f
	f77 -c $(flags) em_2d_slp.f
slp_arc_host.o: slp_arc_host.f
	f77 -c $(flags) slp_arc_host.f
slp_arc_blended.o: slp_arc_blended.f
	f77 -c $(flags) slp_arc_blended.f
em_2d_dlp.o: em_2d_dlp.f
	f77 -c $(flags) em_2d_dlp.f
taylor.o: taylor.f
	f77 -c $(flags) taylor.f
em_2d_rhe.o: em_2d_rhe.f
	f77 -c $(flags) em_2d_rhe.f
em_2d_pdfl.o: em_2d_pdfl.f
	f77 -c $(flags) em_2d_pdfl.f
sgf_2d_fs.o: sgf_2d_fs.f
	f77 -c $(flags) sgf_2d_fs.f
sgf_2d_w.o: sgf_2d_w.f
	f77 -c $(flags) sgf_2d_w.f
sgf_2d_1p_w.o: sgf_2d_1p_w.f
	f77 -c $(flags) sgf_2d_1p_w.f
sgf_2d_1p_ww.o: sgf_2d_1p_ww.f
	f77 -c $(flags) sgf_2d_1p_ww.f
sgf_2d_1p_ww_fga.o: sgf_2d_1p_ww_fga.f
	f77 -c $(flags) sgf_2d_1p_ww_fga.f
sgf_2d_2p.o: sgf_2d_2p.f
	f77 -c $(flags) sgf_2d_2p.f
sgf_2d_2p_ewald.o: sgf_2d_2p_ewald.f
	f77 -c $(flags) sgf_2d_2p_ewald.f
sgf_2d_2p_qqq.o: sgf_2d_2p_qqq.f
	f77 -c $(flags) sgf_2d_2p_qqq.f
sgf_2d_2p_vvv.o: sgf_2d_2p_vvv.f
	f77 -c $(flags) sgf_2d_2p_vvv.f
gauss_leg.o: gauss_leg.f
	f77 -c $(flags) gauss_leg.f
sgf_2d_2p_int.o: sgf_2d_2p_int.f
	f77 -c $(flags) sgf_2d_2p_int.f
prd_2d.o: prd_2d.f
	f77 -c $(flags) prd_2d.f
prd_2d_interp.o: prd_2d_interp.f
	f77 -c $(flags) prd_2d_interp.f
gel_large.o: gel_large.f
	f77 -c $(flags) gel_large.f
gel_small.o: gel_small.f
	f77 -c $(flags) gel_small.f
exp_int.o: exp_int.f
	f77 -c $(flags) exp_int.f
#
# all
# -----
#
all:
	make em_2d
#
# clean
# -----
#
clean:
	rm -f core
	rm -f em_2d $(OBJ)
	rm -f em_2d.rst em_2d.inp em_2d.rhe em_2d.diag em_2d.xy
