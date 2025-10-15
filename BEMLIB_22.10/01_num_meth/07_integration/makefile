#
# filon
# -----
#
filon: filon_dr.o filon.o 
	f77 -o filon filon_dr.o filon.o
filon.o: filon.f
	f77 -c filon.f
#
# gauss_her 
# ---------
#
gauss_her: gauss_her_dr.o gauss_her.o gauss_her_fun.o
	f77 -o gauss_her gauss_her_dr.o gauss_her.o gauss_her_fun.o
gauss_her_dr.o: gauss_her_dr.f
	f77 -c gauss_her_dr.f
gauss_her_fun.o: gauss_her_fun.f
	f77 -c gauss_her_fun.f
gauss_her.o: gauss_her.f
	f77 -c gauss_her.f
#
# gauss_lag
# ---------
#
gauss_lag: gauss_lag_dr.o gauss_lag.o gauss_lag_fun.o
	f77 -o gauss_lag gauss_lag_dr.o gauss_lag.o gauss_lag_fun.o
gauss_lag_dr.o: gauss_lag_dr.f
	f77 -c gauss_lag_dr.f
gauss_lag_fun.o: gauss_lag_fun.f
	f77 -c gauss_lag_fun.f
gauss_lag.o: gauss_lag.f
	f77 -c gauss_lag.f
#
# gauss_log
# ---------
#
gauss_log: gauss_log_dr.o gauss_log.o gauss_log_fun.o
	f77 -o gauss_log gauss_log_dr.o gauss_log.o gauss_log_fun.o
gauss_log_dr.o: gauss_log_dr.f
	f77 -c gauss_log_dr.f
gauss_log_fun.o: gauss_log_fun.f
	f77 -c gauss_log_fun.f
gauss_log.o: gauss_log.f
	f77 -c gauss_log.f
# 
# gauss_sph
# ---------
#
gauss_sph: gauss_sph.f
	f77 -c gauss_sph.f
# 
# gauss_sph_cc
# ------------
#
gauss_sph_cc: gauss_sph_cc_dr.o gauss_sph_cc.o
	g++ -o gauss_sph_cc gauss_sph_cc_dr.o gauss_sph_cc.o
gauss_sph_cc_dr.o: gauss_sph_cc_dr.cc
	g++ -c gauss_sph_cc_dr.cc
gauss_sph_cc.o: gauss_sph_cc.cc
	g++ -c gauss_sph_cc.cc
#
# gauss_trgl
# ---------
#
gauss_trgl: gauss_trgl_dr.o gauss_trgl.o gauss_trgl_fun.o
	f77 -o gauss_trgl gauss_trgl_dr.o gauss_trgl.o gauss_trgl_fun.o
gauss_trgl_dr.o: gauss_trgl_dr.f
	f77 -c gauss_trgl_dr.f
gauss_trgl_fun.o: gauss_trgl_fun.f
	f77 -c gauss_trgl_fun.f
gauss_trgl.o: gauss_trgl.f
	f77 -c gauss_trgl.f
#
# integral_1d
# --------
#
integral_1d: integral_1d.o integral_1d_fun.o gauss_leg.o
	f77 -o integral_1d integral_1d.o integral_1d_fun.o gauss_leg.o
integral_1d.o: integral_1d.f
	f77 -c integral_1d.f
integral_1d_fun.o: integral_1d_fun.f
	f77 -c integral_1d_fun.f
gauss_leg.o: gauss_leg.f
	f77 -c gauss_leg.f
#
# integral_quadri
# ---------------
#
integral_quadri: integral_quadri.o integral_quadri_fun.o 
	f77 -o integral_quadri integral_quadri.o integral_quadri_fun.o
integral_quadri.o: integral_quadri.f
	f77 -c integral_quadri.f
integral_quadri_fun.o: integral_quadri_fun.f
	f77 -c integral_quadri_fun.f
#
# integral_rec
# ------------
#
integral_rec: integral_rec.o integral_rec_fun.o 
	f77 -o integral_rec integral_rec.o integral_rec_fun.o
integral_rec.o: integral_rec.f
	f77 -c integral_rec.f
integral_rec_fun.o: integral_rec_fun.f
	f77 -c integral_rec_fun.f
#
# lobatto
# -------
#
lobatto: lobatto_dr.o lobatto_dr_fun.o lobatto.o
	f77 -o lobatto lobatto_dr.o lobatto_dr_fun.o lobatto.o
lobatto_dr.o: lobatto_dr.f
	f77 -c lobatto_dr.f
lobatto_dr_fun.o: lobatto_dr_fun.f
	f77 -c lobatto_dr_fun.f
gauss_lob_fun.o: gauss_lob_fun.f
	f77 -c gauss_lob_fun.f
lobatto.o: lobatto.f
	f77 -c lobatto.f
#
# radau
# -----
#
radau: radau_dr.o radau_fun.o radau.o
	f77 -o radau radau_dr.o radau_fun.o radau.o
radau_dr.o: radau_dr.f
	f77 -c radau_dr.f
radau_fun.o: radau_fun.f
	f77 -c radau_fun.f
radau.o: radau.f
	f77 -c radau.f
#
# all
# ---
#
all:
	make filon
	make gauss_her
	make gauss_lag
	make gauss_log
	make gauss_sph
	make gauss_sph_cc
	make gauss_trgl
	make integral_1d
	make integral_quadri
	make integral_rec
	make lobatto
	make radau
#
# all_cc
# ------
#
all_cc:
	make gauss_sph_cc
#
# clean
# -----
#
clean:
	rm -f core
	rm -f gauss_her_dr.o gauss_her_fun.o gauss_her.o gauss_her gauss_her.out
	rm -f gauss_lag_dr.o gauss_lag_fun.o gauss_lag.o gauss_lag gauss_lag.out
	rm -f gauss_log_dr.o gauss_log_fun.o gauss_log.o gauss_log gauss_log.out
	rm -f gauss_sph.o 
	rm -f gauss_sph_cc.o gauss_sph_cc_dr.o gauss_sph_cc
	rm -f gauss_trgl_dr.o gauss_trgl_fun.o gauss_trgl.o gauss_trgl gauss_trgl.out
	rm -f lobatto_dr.o lobatto_dr_fun.o lobatto.o lobatto lobatto.out
	rm -f radau_dr.o radau_fun.o radau.o radau radau.out
	rm -f filon_dr.o filon.o filon
	rm -f integral_1d.o integral_1d_fun.o gauss_leg.o integral_1d integral_1d.out
	rm -f integral_rec.o integral_rec_fun.o integral_rec integral_rec.out
	rm -f integral_quadri.o integral_quadri_fun.o integral_quadri integral_quadri.out
