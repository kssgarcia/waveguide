# OBJECTS
#
OBJ1 = prtcl_2d.o prtcl_2d_aux.o prtcl_2d_strml.o prtcl_2d_vel.o 
OBJ2 = prtcl_2d_geo.o prtcl_2d_sdlp.o 
OBJ3 = sgf_2d_fs.o sgf_2d_w.o sgf_2d_1p_w.o sgf_2d_1p_ww.o
OBJ3A= sgf_2d_2p.o sgf_2d_2p_qqq.o sgf_2d_2p_vvv.o sgf_2d_2p_ewald.o 
OBJ4 = exp_int.o gauss_leg.o
OBJ5 = gel_inv.o
OBJ6 = prtcl_2d_sys.o
OBJ7 = contact.o link.o prtcl_2d_sys1.o
OBJs = $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ3A) $(OBJ4) $(OBJ5) $(OBJ6) $(OBJ7)
#
# link
# ----
#
prtcl_2d: $(OBJs)
	f77 -o prtcl_2d $(OBJs)
#
# compile
# ------
#
gauss_leg.o: gauss_leg.f
	f77 -c gauss_leg.f
gel_inv.o: gel_inv.f
	f77 -c gel_inv.f
prtcl_2d.o: prtcl_2d.f
	f77 -c prtcl_2d.f
prtcl_2d_geo.o: prtcl_2d_geo.f
	f77 -c prtcl_2d_geo.f
prtcl_2d_aux.o: prtcl_2d_aux.f
	f77 -c prtcl_2d_aux.f
prtcl_2d_sys.o: prtcl_2d_sys.f
	f77 -c prtcl_2d_sys.f
prtcl_2d_sdlp.o: prtcl_2d_sdlp.f
	f77 -c prtcl_2d_sdlp.f
prtcl_2d_strml.o: prtcl_2d_strml.f
	f77 -c prtcl_2d_strml.f
prtcl_2d_vel.o: prtcl_2d_vel.f
	f77 -c prtcl_2d_vel.f
sgf_2d_fs.o: sgf_2d_fs.f
	f77 -c sgf_2d_fs.f
sgf_2d_w.o: sgf_2d_w.f
	f77 -c sgf_2d_w.f
sgf_2d_1p_w.o: sgf_2d_1p_w.f
	f77 -c sgf_2d_1p_w.f
sgf_2d_2p.o: sgf_2d_2p.f
	f77 -c sgf_2d_2p.f
sgf_2d_2p_qqq.o: sgf_2d_2p_qqq.f
	f77 -c sgf_2d_2p_qqq.f
sgf_2d_2p_vvv.o: sgf_2d_2p_vvv.f
	f77 -c sgf_2d_2p_vvv.f
sgf_2d_2p_ewald.o: sgf_2d_2p_ewald.f
	f77 -c sgf_2d_2p_ewald.f
sgf_2d_1p_ww.o: sgf_2d_1p_ww.f
	f77 -c sgf_2d_1p_ww.f
exp_int.o: exp_int.f
	f77 -c exp_int.f
#
prtcl_2d_sys1.o: prtcl_2d_sys1.f
	f77 -c prtcl_2d_sys1.f
link.o: link.f
	f77 -c link.f
contact.o: contact.f
	f77 -c contact.f
#
#  prtcl_2d_gen generates random distributions for input
# ------------------------------------------------------
#
prtcl_2d_gen: prtcl_2d_gen.o random.o
	f77 -o prtcl_2d_gen prtcl_2d_gen.o random.o
prtcl_2d_gen.o: prtcl_2d_gen.f
	f77 -c prtcl_2d_gen.f
random.o: random.f
	f77 -c random.f
#
# all:
# ---
#
all:
	make prtcl_2d
	make prtcl_2d_gen
#
# clean:
# -----
#
clean:
	rm -f PLOTDAT AAA-xmgr
	rm -f core $(OBJs) prtcl_2d
	rm -f prtcl_2d.inp prtcl_2d.str prtcl_2d.error prtcl_2d.xy
	rm -f prtcl_2d.strml
	rm -f prtcl_2d_gen.o random.o prtcl_2d_gen
#
# purge
# -----
#
purge:
	rm -f PLOTDAT AAA-xmgr
	rm -f core $(OBJ)
	rm -f prtcl_2d.inp prtcl_2d.str prtcl_2d.error prtcl_2d.xy
	rm -f prtcl_2d.strml
	rm -f prtcl_2d_gen.o random.o 
#
# tar
# ---
#
tar:
	tar cvf prtcl_2d.tar makefile *.f *.dat
