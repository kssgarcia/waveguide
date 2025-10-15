#
# chol
# ---
#
chol: chol_dr.o chol_c.o chol_r.o 
	f77 -o chol chol_dr.o  chol_c.o chol_r.o 
chol_dr.o: chol_dr.f
	f77 -c chol_dr.f
chol_c.o: chol_c.f
	f77 -c chol_c.f
chol_r.o: chol_r.f
	f77 -c chol_r.f
#
# det_hess
# --------
#
det_hess: det_hess.o det_hess_l.o det_hess_u.o 
	f77 -o det_hess det_hess.o det_hess_l.o det_hess_u.o 
det_hess.o: det_hess.f
	f77 -c det_hess.f
det_hess_l.o: det_hess_l.f
	f77 -c det_hess_l.f
det_hess_u.o: det_hess_u.f
	f77 -c det_hess_u.f
#
# gram_schmidt.f
# -------
#
gram_schmidt: gram_schmidt.o 
	f77 -o gram_schmidt gram_schmidt.o 
#
# inv_l
# -----
#
inv_l: inv_l_dr.o inv_l.o
	f77 -o inv_l inv_l_dr.o inv_l.o 
inv_l_dr.o: inv_l_dr.f
	f77 -c inv_l_dr.f
inv_l.o: inv_l.f
	f77 -c inv_l.f
#
# inv_u
# -----
#
inv_u: inv_u_dr.o inv_u.o
	f77 -o inv_u inv_u_dr.o inv_u.o 
inv_u_dr.o: inv_u_dr.f
	f77 -c inv_u_dr.f
inv_u.o: inv_u.f
	f77 -c inv_u.f
#
# ldu
# ---
#
ldu: ldu_dr.o ldu.o
	f77 -o ldu ldu_dr.o ldu.o 
ldu_dr.o: ldu_dr.f 
	f77 -c ldu_dr.f 
ldu.o: ldu.f 
	f77 -c ldu.f 
#
# lu_c
# ----
#
lu_c: lu_c_dr.o lu_c.o inv_l.o inv_u.o
	f77 -o lu_c lu_c_dr.o lu_c.o inv_l.o inv_u.o
lu_c_dr.o: lu_c_dr.f
	f77 -c lu_c_dr.f
lu_c.o: lu_c.f
	f77 -c lu_c.f
#
# lu_d
# ----
#
lu_d: lu_d_dr.o lu_d.o inv_l.o inv_u.o
	f77 -o lu_d lu_d_dr.o lu_d.o inv_l.o inv_u.o
lu_d_dr.o: lu_d_dr.f
	f77 -c lu_d_dr.f
lu_d.o: lu_d.f
	f77 -c lu_d.f
#
# lu_d_t
# ------
#
lu_d_t: lu_d_t_dr.o lu_d_t.o 
	f77 -o lu_d_t lu_d_t_dr.o lu_d_t.o 
lu_d_t_dr.o: lu_d_t_dr.f
	f77 -c lu_d_t_dr.f
lu_d_t.o: lu_d_t.f
	f77 -c lu_d_t.f
#
# matrix_mult
# -----------
#
matrix_mult: matrix_mult.o
	f77 -o matrix_mult matrix_mult.f
#
# matrix_pow
# -----------
#
matrix_pow: matrix_pow.o
	f77 -o matrix_pow matrix_pow.f
#
# qr_dec
# ------
#
qr_dec: qr_dec.o qr_gs.o qr_rot.o qr_ref.o
	f77 -o qr_dec qr_dec.o qr_gs.o qr_rot.o qr_ref.o
qr_dec.o.o: qr_dec.f
	f77 -c qr_dec.f
qr_gs.o: qr_gs.f
	f77 -c qr_gs.f
qr_rot.o: qr_rot.f
	f77 -c qr_rot.f
qr_ref.o: qr_ref.f
	f77 -c qr_ref.f
#
# all
# ---
#
all:
	make chol
	make det_hess
	make gram_schmidt
	make inv_u
	make ldu
	make lu_c
	make lu_d
	make lu_d_t
	make matrix_mult
	make matrix_pow
	make qr_dec
#
# clean
# -----
#
clean:
	rm -f chol_dr.o chol_c.o chol_r.o chol 
	rm -f det_hess.o det_hess_l.o det_hess_u.o det_hess
	rm -f gram_sc.o gram_sc
	rm -f inv_l_dr.o inv_l.o inv_l
	rm -f inv_u_dr.o inv_u.o inv_u
	rm -f ldu_dr.o ldu.o ldu
	rm -f lu_c_dr.o lu_c.o inv_l.o inv_u.o lu_c
	rm -f lu_d_dr.o lu_d.o inv_l.o inv_u.o lu_d
	rm -f lu_d_t_dr.o lu_d_t.o lu_d_t
	rm -f matrix_mult.o matrix_mult
	rm -f matrix_pow.o matrix_pow
	rm -f qr_dec.o qr_gs.o qr_rot.o qr_ref.o qr_dec
