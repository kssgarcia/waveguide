OBJg = 
#
OBJ1 = chan2l.o chan2l_slv.o gauss_leg.o chan2l_sdlp.o 
OBJ2 = gel.o thomas_pr.o chan2l_cd.o
OBJ3 = sgf_2d_fs.o sgf_2d_1p_ww.o sgf_2d_1p_ww_fga.o
OBJ4 = splc_mnmx.o splc_pr.o splc_pr_geo.o
OBJ5 = prd_2d_pr_splc.o prd_2d_pr_splc_interp.o
OBJ  = $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5)
#
#
# libraries
# ---------
#
LIBS = -lm 
#
# compiler flags
# --------------
#
FLAGS =
#
# link
# ----
#
chan2l: $(OBJ)
	f77 -o chan2l $(OBJ) $(LIBS)
#
# compile
# -------
#
chan2l.o: chan2l.f
	f77 -c chan2l.f
chan2l_sdlp.o: chan2l_sdlp.f
	f77 -c chan2l_sdlp.f
chan2l_slv.o: chan2l_slv.f
	f77 -c chan2l_slv.f
gauss_leg.o: gauss_leg.f
	f77 -c gauss_leg.f
redist.o: redist.f
	f77 -c redist.f
gel.o: gel.f
	f77 -c gel.f
prd_2d_pr_splc.o: prd_2d_pr_splc.f
	f77 -c prd_2d_pr_splc.f
prd_2d_pr_splc_interp.o: prd_2d_pr_splc_interp.f
	f77 -c prd_2d_pr_splc_interp.f
chan2l_cd.o: chan2l_cd.f
	f77 -c chan2l_cd.f
splc_pr.o: splc_pr.f
	f77 -c splc_pr.f
splc_mnmx.o: splc_mnmx.f
	f77 -c splc_mnmx.f
splc_pr_geo.o: splc_pr_geo.f
	f77 -c splc_pr_geo.f
sgf_2d_fs.o: sgf_2d_fs.f
	f77 -c sgf_2d_fs.f
sgf_2d_1p_ww.o: sgf_2d_1p_ww.f
	f77 -c sgf_2d_1p_ww.f
sgf_2d_1p_ww_fga.o: sgf_2d_1p_ww_fga.f
	f77 -c sgf_2d_1p_ww_fga.f
thomas_pr.o: thomas_pr.f
	f77 -c $(FLAGS) thomas_pr.f
#
# compile graphics
# ----------------
#
#
# all
# ---
#
all:
	make chan2l
#
# clean
# -----
#
clean:
	rm -f core $(OBJ)
	rm -f chan2l.xy chan2l.diag chan2l.rst
	rm -f chan2l
#
# purge
# -----
#
purge:
	rm -f core $(OBJ)
	rm -f chan2l.xy chan2l.diag chan2l.rst
