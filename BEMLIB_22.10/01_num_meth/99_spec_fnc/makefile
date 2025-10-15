#
# special functions
# =================
#
OBJ1 = spec_fnc.o ell_int.o erfun.o erfuni.o 
OBJ2 = bess_I01K01.o bess_J0.o bess_J1.o ber_bei_0.o
OBJ3 = ker_kei_0.o ker_kei_0_large.o
OBJ  = $(OBJ1) $(OBJ2) $(OBJ3)
#
spec_fnc: $(OBJ)
	f77 -o spec_fnc $(OBJ)
spec_fnc.o: spec_fnc.f
	f77 -c spec_fnc.f
ell_int.o: ell_int.f
	f77 -c ell_int.f
erfun.o: erfun.f
	f77 -c erfun.f
erfuni.o: erfuni.f
	f77 -c erfuni.f
bess_J0.o: bess_J0.f
	f77 -c bess_J0.f
bess_J1.o: bess_J1.f
	f77 -c bess_J1.f
ber_bei_0.o: ber_bei_0.f
	f77 -c ber_bei_0.f
ker_kei_0.o: ker_kei_0.f
	f77 -c ker_kei_0.f
ker_kei_0_large.o: ker_kei_0_large.f
	f77 -c ker_kei_0_large.f
bess_I01K01.o: bess_I01K01.f
	f77 -c bess_I01K01.f
#
# poly_ortho
# ==========
#
OBJ001 = poly_ortho.o poly_ortho_her.o poly_ortho_lgr.o poly_ortho_lin.o
OBJ002 = poly_ortho_log.o poly_ortho_rad.o poly_ortho_leg.o
OBJ003 = poly_ortho_G.o poly_ortho_prj.o
OBJ000 = $(OBJ001) $(OBJ002) $(OBJ003)
#
poly_ortho: $(OBJ000)
	f77 -o poly_ortho $(OBJ000)
poly_ortho_leg: poly_ortho_leg.f
	f77 -c poly_ortho_leg.f
poly_ortho_rad: poly_ortho_rad.f
	f77 -c poly_ortho_rad.f
poly_ortho_her: poly_ortho_her.f
	f77 -c poly_ortho_her.f
poly_ortho_lgr: poly_ortho_lgr.f
	f77 -c poly_ortho_lgr.f
poly_ortho_lin: poly_ortho_lin.f
	f77 -c poly_ortho_lin.f
poly_ortho_log: poly_ortho_log.f
	f77 -c poly_ortho_log.f
poly_ortho_G: poly_ortho_G.f
	f77 -c poly_ortho_G.f
#
#
#
leg_poly: leg_poly.o leg_poly_dr.o
	f77 -o leg_poly leg_poly.o leg_poly_dr.o
leg_poly_dr.o: leg_poly_dr.f
	f77 -c leg_poly_dr.f
leg_poly.o: leg_poly.f
	f77 -c leg_poly.f
#
# compile c++
# -----------
#
ell_int_cc: ell_int_cc.o
	g++ -c ell_int_cc.cc
exp_int_cc: exp_int_cc.o
	g++ -c exp_int_cc.cc
#
#
# all
# ---
#
all:
	make spec_fnc
	make poly_ortho
	make ell_int_cc
	make exp_int_cc
#
# clean
# -----
#
clean:
	rm -f core spec_fnc.out
	rm -f $(OBJ) poly_ortho.o poly_ortho
	rm -f spec_fnc
	rm -f ell_int.o exp_int.o
	rm -f poly_ortho $(OBJ000)
	rm -f leg_poly leg_poly.o leg_poly_dr.o
#
#
# purge
# --------
#
purge:
	rm -f core
	rm -f $(OBJ) poly_ortho.o
	rm -f ell_int.o exp_int.o
	rm -f ell_int_cc.o exp_int_cc.o
	rm -f $(OBJ000)
	rm -f leg_poly.o leg_poly_dr.o
