#
# object files
#
OBJ1 = caps_3d.o trgl6_octa.o trgl6_icos.o
OBJ11= crvm_3d.o crvm_3d_interp.o
OBJ2 = sgf_3d_fs.o sgf_3d_w.o
OBJ2A= sgf_3d_3p.o sgf_3d_3p_ewald.o sgf_3d_3p_qqq.o sgf_3d_3p_vvv.o
OBJ3 = abc.o interp_p.o xy_slice.o taylor.o inclination.o printel.o
OBJ32 = elm_geom.o
OBJ31= rheology.o rheology_interp.o
OBJ33 = caps_3d_dfel.o caps_3d_dfel_interp.o
OBJ4 = caps_3d_vel.o deflation.o 
OBJ41= caps_3d_slp.o caps_3d_slp_integrate.o caps_3d_slp_sing.o caps_3d_slp_interp.o 
OBJ42= sdlp_3d.o sdlp_3d_interp.o sdlp_3d_integral.o
OBJ43 = sgrad_3d.o sgrad_3d_interp.o 
OBJ5 = elten_3d.o elten_3d_interp.o vtg_3d.o vtg_3d_interp.o
OBJ6 = cramer_33.o gauss_trgl.o gauss_leg.o ttr.o
OBJA  = $(OBJ1) $(OBJ11) $(OBJ2) $(OBJ3) $(OBJ31) $(OBJ4) $(OBJ33)
OBJB  = $(OBJ41) $(OBJ42) $(OBJ5) $(OBJ6) $(OBJ32) $(OBJ43)
OBJ   = $(OBJA) $(OBJB) $(OBJ2A)
#
# link
#
caps_3d: $(OBJ)
	f77 -o caps_3d $(OBJ)
#
# compile
#
caps_3d.o: caps_3d.f
	f77 -c caps_3d.f 
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
printel.o: printel.f 
	f77 -c printel.f 
abc.o: abc.f
	f77 -c abc.f
interp_p.o: interp_p.f
	f77 -c interp_p.f
crvm_3d.o: crvm_3d.f 
	f77 -c crvm_3d.f 
crvm_3d_interp.o: crvm_3d_interp.f 
	f77 -c crvm_3d_interp.f 
elm_geom.o: elm_geom.f
	f77 -c elm_geom.f
caps_3d_dfel.o: caps_3d_dfel.f 
	f77 -c caps_3d_dfel.f 
caps_3d_dfel_interp.o: caps_3d_dfel_interp.f 
	f77 -c caps_3d_dfel_interp.f 
rheology.o: rheology.f 
	f77 -c rheology.f 
rheology_interp.o: rheology_interp.f 
	f77 -c  rheology_interp.f 
xy_slice.o: xy_slice.f
	f77 -c xy_slice.f
inclination.o: inclination.f
	f77 -c inclination.f
taylor.o: taylor.f
	f77 -c taylor.f
caps_3d_vel.o: caps_3d_vel.f 
	f77 -c caps_3d_vel.f 
deflation.o: deflation.f
	f77 -c deflation.f 
caps_3d_slp.o: caps_3d_slp.f 
	f77 -c caps_3d_slp.f 
caps_3d_slp_integrate.o: caps_3d_slp_integrate.f 
	f77 -c caps_3d_slp_integrate.f 
caps_3d_slp_interp.o: caps_3d_slp_interp.f 
	f77 -c caps_3d_slp_interp.f 
caps_3d_slp_sing.o: caps_3d_slp_sing.f 
	f77 -c caps_3d_slp_sing.f 
sdlp_3d.o: sdlp_3d.f
	f77 -c sdlp_3d.f 
sdlp_3d_interp.o: sdlp_3d_interp.f
	f77 -c sdlp_3d_interp.f
sdlp_3d_integral.o: sdlp_3d_integral.f
	f77 -c sdlp_3d_integral.f
elten_3d.o: elten_3d.f
	f77 -c elten_3d.f
elten_3d_interp.o: elten_3d_interp.f
	f77 -c elten_3d_interp.f
vtg_3d.o: vtg_3d.f
	f77 -c vtg_3d.f
vtg_3d_interp.o: vtg_3d_interp.f
	f77 -c vtg_3d_interp.f
sgrad_3d.o: sgrad_3d.f
	f77 -c sgrad_3d.f
sgrad_3d_interp.o: sgrad_3d_interp.f
	f77 -c sgrad_3d_interp.f
cramer_33.o: cramer_33.f
	f77 -c cramer_33.f
gauss_leg.o: gauss_leg.f 
	f77 -c gauss_leg.f 
gauss_trgl.o: gauss_trgl.f 
	f77 -c gauss_trgl.f 
ttr.o: ttr.f
	f77 -c ttr.f
#
# all
# ---
#
all:
	make caps_3d
#
# tar
# ---
#
tar:
	tar cvf caps_3d.tar *.f caps_3d.dat makefile
#
# clean
# -----
#
clean:
	rm -f core PLOTDAT AAA-xmgr
	rm -f $(OBJ) caps_3d
	rm -f caps_3d.xy caps_3d.xyz caps_3d.diag
	rm -f caps_3d.rst caps_3d.inp
