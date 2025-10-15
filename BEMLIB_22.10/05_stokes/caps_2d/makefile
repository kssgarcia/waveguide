#
# object files
#
OBJ1 = caps_2d.o caps_2d_vel.o caps_2d_dlp.o
OBJ2 = arc_2d.o gel.o
OBJ3 = sgf_2d_fs.o sgf_2d_w.o sgf_2d_1p_w.o
OBJ3A = sgf_2d_1p_ww.o sgf_2d_1p_ww_fga.o
OBJ3B =  sgf_2d_2p.o sgf_2d_2p_ewald.o sgf_2d_2p_qqq.o sgf_2d_2p_vvv.o
OBJ4 = sgf_2d_2p_int.o gauss_leg.o exp_int.o
OBJ5 = taylor.o crv_max.o deflation.o caps_2d_rhe.o caps_2d_pdfl.o 
OBJ6 = verbal.o expand.o
OBJ7 = prd_2d.o prd_2d_interp.o prd_2d_qs.o prd_2d_qs_interp.o
OBJ8 = caps_2d_slp.o slp_arc_host.o slp_arc_blended.o
OBJ9 = caps_2d_cd.o interp_en.o interp_mn.o
OBJA = $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ8)
OBJB = $(OBJ6) $(OBJ7) $(OBJ3A) $(OBJ3B) $(OBJ9)
#
OBJ  = $(OBJA) $(OBJB)
#
# libraries
# ---------
#
LIBS =
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
caps_2d: $(OBJ)
	f77 -o caps_2d $(OBJ) $(LIBS)
#
# compile
# -------
#
caps_2d.o: caps_2d.f
	f77 -c $(flags) caps_2d.f
arc_2d.o: arc_2d.f
	f77 -c $(flags) arc_2d.f
caps_2d_cd.o: caps_2d_cd.f
	f77 -c $(flags) caps_2d_cd.f
interp_en.o: interp_en.f
	f77 -c $(flags) interp_en.f
interp_mn.o: interp_mn.f
	f77 -c $(flags) interp_mn.f
caps_2d_vel.o: caps_2d_vel.f
	f77 -c $(flags) caps_2d_vel.f
caps_2d_slp.o: caps_2d_slp.f
	f77 -c $(flags) caps_2d_slp.f
slp_arc_host.o: slp_arc_host.f
	f77 -c $(flags) slp_arc_host.f
slp_arc_blended.o: slp_arc_blended.f
	f77 -c $(flags) slp_arc_blended.f
caps_2d_dlp.o: caps_2d_dlp.f
	f77 -c $(flags) caps_2d_dlp.f
taylor.o: taylor.f
	f77 -c $(flags) taylor.f
crv_max.o: crv_max.f
	f77 -c $(flags) crv_max.f
deflation.o: deflation.f
	f77 -c $(flags) deflation.f
caps_2d_rhe.o: caps_2d_rhe.f
	f77 -c $(flags) caps_2d_rhe.f
caps_2d_pdfl.o: caps_2d_pdfl.f
	f77 -c $(flags) caps_2d_pdfl.f
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
sgf_2d_2p_int.o: sgf_2d_2p_int.f
	f77 -c $(flags) sgf_2d_2p_int.f
gauss_leg.o: gauss_leg.f
	f77 -c $(flags) gauss_leg.f
prd_2d.o: prd_2d.f
	f77 -c $(flags) prd_2d.f
prd_2d_interp.o: prd_2d_interp.f
	f77 -c $(flags) prd_2d_interp.f
prd_2d_qs.o: prd_2d_qs.f
	f77 -c $(flags) prd_2d_qs.f
prd_2d_qs_interp.o: prd_2d_qs_interp.f
	f77 -c $(flags) prd_2d_qs_interp.f
gel.o: gel.f
	f77 -c $(flags) gel.f
exp_int.o: exp_int.f
	f77 -c $(flags) exp_int.f
caps_2d_gr.o: caps_2d_gr.f
	f77 -c $(flags) caps_2d_gr.f
verbal.o: verbal.f
	f77 -c $(flags) verbal.f
expand.o: expand.f
	f77 -c $(flags) expand.f
#
# all
# ---
#
all:
	make caps_2d
#
# clean
# -----
#
clean:
	rm -f $(OBJ) core
	rm -f $(OBJ) caps_2d
	rm -f caps_2d.diag caps_2d.xy caps_2d.rhe caps_2d.rst
#
# tar
# -----
#
tar:
	tar cvf caps_2d.tar *.f *.dat makefile
