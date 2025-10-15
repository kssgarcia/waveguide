#
# bcg
# ---
#
#
bcg: bcg_dr.o bcg.o
	f77 -o bcg bcg_dr.o bcg.o
bcg_dr.o: bcg_dr.f
	f77 -c bcg_dr.f
bcg.o: bcg.f
	f77 -c bcg.f
#
# cg
# --
#
cg: cg_dr.o cg.o
	f77 -o cg cg_dr.o cg.o
cg_dr.o: cg_dr.f
	f77 -c cg_dr.f
cg.o: cg.f
	f77 -c cg.f
#
# cycl_red
# --------
#
cycl_red: cycl_red_dr.o cycl_red.o
	f77 -o cycl_red cycl_red_dr.o cycl_red.o
cycl_red_dr.o: cycl_red_dr.f
	f77 -c cycl_red_dr.f
cycl_red.o: cycl_red.f
	f77 -c cycl_red.f
#
# gel
# ---
#
gel: gel_dr.o gel.o
	f77 -o gel gel_dr.o gel.o
gel_dr.o: gel_dr.f
	f77 -c gel_dr.f
gel.o: gel.f
	f77 -c gel.f
#
# gel_inv
# -------
#
gel_inv: gel_inv_dr.o gel_inv.o
	f77 -o gel_inv gel_inv_dr.o gel_inv.o
gel_inv_dr.o: gel_inv_dr.f
	f77 -c gel_inv_dr.f
gel_inv.o: gel_inv.f
	f77 -c gel_inv.f
#
# gel_mrhs
# ------
#
gel_mrhs: gel_mrhs_dr.o gel_mrhs.o
	f77 -o gel_mrhs gel_mrhs_dr.o gel_mrhs.o
gel_mrhs_dr.o:
	f77 -c gel_mrhs_dr.f
gel_mrhs.o:
	f77 -c gel_mrhs.f
#
# thomas
# ------
#
thomas: thomas_dr.o thomas.o thomas_c.o thomas_pr.o
	f77 -o thomas thomas_dr.o thomas.o thomas_c.o thomas_pr.o
thomas_dr.o: thomas_dr.f
	f77 -c thomas_dr.f
thomas.o: thomas.f
	f77 -c thomas.f
thomas_c.o: thomas_c.f
	f77 -c thomas_c.f
thomas_pr.o: thomas_pr.f
	f77 -c thomas_pr.f
#
# penta
# -----
#
penta: penta_dr.o penta.o
	f77 -o penta penta_dr.o penta.o
penta_dr.o:
	f77 -c penta_dr.f
penta.o:
	f77 -c penta.f
#
# pois_fds_DNDD
# --------------
#
pois_fds_DNDD: pois_fds_DNDD_dr.o pois_fds_DNDD.o
	f77 -o pois_fds_DNDD pois_fds_DNDD_dr.o pois_fds_DNDD.o
pois_fds_DNDD_dr.o: pois_fds_DNDD_dr.f  
	f77 -c pois_fds_DNDD_dr.f
pois_fds_DNDD.o: pois_fds_DNDD.f
	f77 -c pois_fds_DNDD.f
#
# all
# ---
#
all:
	make bcg
	make cg
	make cycl_red
	make gel
	make gel_inv
	make gel_mrhs
	make thomas
	make penta
	make pois_fds_DNDD
#
# clean
# -----
#
clean:
	rm -f bcg.o bcg_dr.o bcg
	rm -f cg.o cg_dr.o cg
	rm -f cycl_red_dr.o cycl_red.o cycl_red
	rm -f gel.o gel_dr.o gel
	rm -f gel_inv.o gel_inv_dr.o gel_inv
	rm -f gel_mrhs.o gel_mrhs_dr.o gel_mrhs
	rm -f thomas_dr.o thomas.o thomas_c.o thomas_pr.o thomas
	rm -f penta.o penta_dr.o penta
	rm -f pois_fds_DNDD_dr.o pois_fds_DNDD.o pois_fds_DNDD
#
# purge
# -----
#
purge:
	rm -f bcg.o bcg_dr.o
	rm -f cg.o cg_dr.o
	rm -f cycl_red_dr.o cycl_red.o
	rm -f gel.o gel_dr.o 
	rm -f gel_inv.o gel_inv_dr.o
	rm -f gel_mrhs.o gel_mrhs_dr.o
	rm -f thomas_dr.o thomas.o thomas_c.o thomas_pr.o 
	rm -f penta.o penta_dr.o
