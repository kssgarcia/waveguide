#
# object files
# ------------
#
OBJ1 = drop_3d.o verbal.o trgl6_octa.o trgl6_icos.o
OBJ11= crvm_3d.o crvm_3d_interp.o
OBJ2 = sgf_3d_fs.o sgf_3d_w.o
OBJ2A= sgf_3d_3p.o sgf_3d_3p_ewald.o sgf_3d_3p_qqq.o sgf_3d_3p_vvv.o
OBJ3 = abc.o inclination.o taylor.o xy_slice.o interp_p.o
OBJ39= elm_geom.o printel.o
OBJ33= rheology.o rheology_interp.o
OBJ4 = vel.o deflation.o
OBJ41= sslp_3d.o sslp_3d_integral.o sslp_3d_integral_sing.o sslp_3d_interp.o
OBJ42= sdlp_3d.o sdlp_3d_interp.o sdlp_3d_integral.o 
OBJ5 = sgrad_3d.o sgrad_3d_interp.o
OBJ6 = cramer_33.o gauss_trgl.o gauss_leg.o
OBJ7 = dilt_3d.o dilt_3d_interp.o gel.o
OBJ71= surfact_fem.o surfact_fem_interp.o
OBJ72= surfact_fvm.o surfact_fvm_interp.o surfact_fvm_interp1.o surfact_fvm_etn.o 
OBJA = $(OBJ1) $(OBJ11) $(OBJ2) $(OBJ2A) $(OBJ3) $(OBJ33) $(OBJ41) $(OBJ42)
OBJB = $(OBJ4) $(OBJ5) $(OBJ6) $(OBJ7) $(OBJ39) $(OBJ71) $(OBJ72)
OBJ = $(OBJA) $(OBJB) 
#
#
#
OBJJ1 = surfact.o trgl6_octa.o trgl6_icos.o
OBJJ11= crvm_3d.o crvm_3d_interp.o
OBJJ2 = dilt_3d.o dilt_3d_interp.o gel.o
OBJJ3 = abc.o xy_slice.o interp_p.o
OBJJ39= elm_geom.o printel.o
OBJJ4 = thomas.o
OBJJ5 = sgrad_3d.o sgrad_3d_interp.o
OBJJ9 = gauss_trgl.o cramer_33.o
OBJJ71= surfact_fem.o surfact_fem_interp.o
OBJJ72= surfact_fvm.o surfact_fvm_interp.o surfact_fvm_interp1.o surfact_fvm_etn.o 
OBJJA= $(OBJJ1) $(OBJJ11) $(OBJJ2) $(OBJJ3) $(OBJJ4) $(OBJJ9) $(OBJJ39)
OBJJB= $(OBJJ5) $(OBJJ71) $(OBJJ72)
OBJJ = $(OBJJA) $(OBJJB) 
#
# link
# ----
#
drop_3d: $(OBJ)
	f77 -o drop_3d $(OBJ)
surfact: $(OBJJ)
	f77 -o surfact $(OBJJ)
#
# compile
# -------
#
cramer_33.o: cramer_33.f
	f77 -c cramer_33.f
crvm_3d.o: crvm_3d.f 
	f77 -c crvm_3d.f 
crvm_3d_interp.o: crvm_3d_interp.f 
	f77 -c crvm_3d_interp.f 
dilt_3d.o: dilt_3d.f 
	f77 -c dilt_3d.f 
dilt_3d_interp.o: dilt_3d_interp.f 
	f77 -c dilt_3d_interp.f 
drop_3d.o: drop_3d.f 
	f77 -c drop_3d.f 
elm_geom.o: elm_geom.f 
	f77 -c elm_geom.f 
abc.o: abc.f
	f77 -c abc.f
printel.o: printel.f
	f77 -c printel.f
interp_p.o: interp_p.f
	f77 -c interp_p.f
taylor.o: taylor.f
	f77 -c taylor.f
inclination.o: inclination.f
	f77 -c inclination.f
xy_slice.o: xy_slice.f
	f77 -c xy_slice.f
rheology.o: rheology.f 
	f77 -c  rheology.f 
rheology_interp.o: rheology_interp.f 
	f77 -c  rheology_interp.f 
vel.o: vel.f 
	f77 -c vel.f 
deflation.o: deflation.f 
	f77 -c deflation.f 
sslp_3d.o: sslp_3d.f 
	f77 -c sslp_3d.f 
sslp_3d_interp.o: sslp_3d_interp.f 
	f77 -c sslp_3d_interp.f 
sslp_3d_integral.o: sslp_3d_integral.f 
	f77 -c sslp_3d_integral.f 
sslp_3d_integral_sing.o: sslp_3d_integral_sing.f 
	f77 -c sslp_3d_integral_sing.f 
sdlp_3d.o: sdlp_3d.f 
	f77 -c sdlp_3d.f 
sdlp_3d_interp.o: sdlp_3d_interp.f
	f77 -c sdlp_3d_interp.f
sdlp_3d_integral.o: sdlp_3d_integral.f
	f77 -c sdlp_3d_integral.f
gauss_leg.o: gauss_leg.f 
	f77 -c gauss_leg.f 
gauss_trgl.o: gauss_trgl.f 
	f77 -c gauss_trgl.f 
gel.o: gel.f
	f77 -c gel.f
trgl6_octa.o: trgl6_octa.f 
	f77 -c trgl6_octa.f 
trgl6_icos.o: trgl6_icos.f 
	f77 -c trgl6_icos.f 
sgf_3d_fs.o: sgf_3d_fs.f 
	f77 -c sgf_3d_fs.f 
sgf_3d_w.o: sgf_3d_w.f 
	f77 -c sgf_3d_w.f 
sgf_3d_3p.o: sgf_3d_3p.f 
	f77 -c sgf_3d_3p.f 
sgf_3d_3p_ewald.o: sgf_3d_3p_ewald.f 
	f77 -c sgf_3d_3p_ewald.f 
sgf_3d_3p_qqq.o: sgf_3d_3p_qqq.f 
	f77 -c sgf_3d_3p_qqq.f 
sgf_3d_3p_vvv.o: sgf_3d_3p_vvv.f 
	f77 -c sgf_3d_3p_vvv.f 
sgrad_3d.o: sgrad_3d.f
	f77 -c sgrad_3d.f
sgrad_3d_interp.o: sgrad_3d_interp.f
	f77 -c sgrad_3d_interp.f
surfact_fem.o: surfact_fem.f
	f77 -c surfact_fem.f
surfact_fem_interp.o: surfact_fem_interp.f
	f77 -c surfact_fem_interp.f
surfact_fvm.o: surfact_fvm.f
	f77 -c surfact_fvm.f
surfact_fvm_etn.o: surfact_fvm_etn.f
	f77 -c surfact_fvm_etn.f
surfact_fvm_interp.o: surfact_fvm_interp.f
	f77 -c surfact_fvm_interp.f
surfact_fvm_interp1.o: surfact_fvm_interp1.f
	f77 -c surfact_fvm_interp1.f
verbal.o: verbal.f 
	f77 -c verbal.f 
thomas.o: thomas.f 
	f77 -c thomas.f 
dfel_3d.o: dfel_3d.f 
	f77 -c dfel_3d.f 
dfel_3d_interp.o: dfel_3d_interp.f 
	f77 -c dfel_3d_interp.f 
#
surfact.o: surfact.f
	f77 -c surfact.f
#
# clean
# -----
#
clean:
	rm -f core PLOTDAT AAA-xmgr
	rm -f $(OBJ) drop_3d
	rm -f drop_3d.xy drop_3d.xyz drop_3d.diag
	rm -f drop_3d.rst
	rm -f drop_3d.rhe
	rm -f drop_3d.net
	rm -f $(OBJJ) surfact
	rm -f surfact.xy surfact.xyz surfact.rst
	rm -f surfact.comp surfact.axis
#
# purge
# -----
#
purge:
	rm -f core PLOTDAT AAA-xmgr
	rm -f $(OBJ)
	rm -f drop_3d.xy drop_3d.xyz drop_3d.diag
#
# all
# ---
#
all:
	make drop_3d
	make surfact
#
# tar
# ---
#
tar:
	tar cvf drop_3d.tar makefile *.f *.dat *.m transfer *.m
