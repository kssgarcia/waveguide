#
# char_poly
# --------
#
char_poly: char_poly.o leverrier.o krylov.o gel_krylov.o
	f77 -o char_poly char_poly.o leverrier.o krylov.o gel_krylov.o
leverrier.o: leverrier.f
	f77 -c leverrier.f
krylov.o:  krylov.f
	f77 -c krylov.f
gel_krylov.o: gel_krylov.f 
	f77 -c gel_krylov.f 
	
#
# cond_num
# --------
#
cond_num: cond_num_dr.o cond_num.o power.o gel_power.o
	f77 -o cond_num cond_num_dr.o cond_num.o power.o gel_power.o
cond_num_dr.o: cond_num_dr.f
	f77 -c cond_num_dr.f
cond_num.o: cond_num.f
	f77 -c cond_num.f
#
# -------
# eigen33
# -------
#
eigen33: eigen33.o cubic.o
	f77 -o eigen33 eigen33.o cubic.o
eigen33.o: eigen33.f
	f77 -c eigen33.f
cubic.o: cubic.f
	f77 -c cubic.f
#
# house_sym
# ---------
#
house_sym: house_sym_dr.o house_sym.o 
	f77 -o house_sym house_sym_dr.o house_sym.o 
house_sym_dr.o: house_sym_dr.f
	f77 -c house_sym_dr.f
house_sym.o: house_sym.f
	f77 -c house_sym.f
#
# jacobi
# ------
#
jacobi: jacobi_dr.o jacobi.o 
	f77 -o jacobi jacobi_dr.o jacobi.o 
jacobi_dr.o: jacobi_dr.f
	f77 -c jacobi_dr.f
jacobi.o: jacobi.f
	f77 -c jacobi.f
#
# lanczos
# -------
#
lanczos: lanczos_dr.o lanczos.o 
	f77 -o lanczos lanczos_dr.o lanczos.o 
lanczos_dr.o: lanczos_dr.f
	f77 -c lanczos_dr.f
lanczos.o: lanczos.f
	f77 -c lanczos.f
#
# lanczos_sym
# -----------
#
lanczos_sym: lanczos_sym_dr.o lanczos_sym.o 
	f77 -o lanczos_sym lanczos_sym_dr.o lanczos_sym.o 
lanczos_sym_dr.o: lanczos_sym_dr.f
	f77 -c lanczos_sym_dr.f
lanczos_sym.o: lanczos_sym.f
	f77 -c lanczos_sym.f
#
# power
# ----
#
power: power_dr.o power.o gel_power.o
	f77 -o power power_dr.o power.o gel_power.o
power_dr.o: power_dr.f
	f77 -c power_dr.f
power.o: power.f
	f77 -c power.f
gel_power.o: gel_power.f 
	f77 -c gel_power.f 
#
# power_sym
# ---------
#
power_sym: power_sym_dr.o power_sym.o 
	f77 -o power_sym power_sym_dr.o power_sym.o 
power_sym_dr.o: power_sym_dr.f
	f77 -c power_sym_dr.f
power_sym.o: power_sym.f
	f77 -c power_sym.f
#
# ------
# qr
# ------
#
qr: qr.o qr_rot.o
	f77 -o qr qr.o qr_rot.o
qr.o.o: qr.f
	f77 -c qr.f
qr_rot.o: qr_rot.f
	f77 -c qr_rot.f
#
# all
# ---
#
all:
	make char_poly
	make cond_num
	make eigen33
	make house_sym
	make jacobi
	make lanczos
	make lanczos_sym
	make power
	make power_sym
	make qr
#
# clean
# -----
#
clean:
	rm -f char_poly.o leverrier.o krylov.o gel_krylov.o char_poly
	rm -f cond_num_dr.o cond_num.o power.o gel_power.o cond_num
	rm -f eigen33.o cubic.o eigen33
	rm -f house_sym_dr.o house_sym.o house_sym
	rm -f jacobi_dr.o jacobi.o jacobi
	rm -f lanczos_dr.o lanczos.o lanczos
	rm -f lanczos_sym_dr.o lanczos_sym.o lanczos_sym
	rm -f power_dr.o power.o gel_power.o power
	rm -f power_sym_dr.o power_sym.o power_sym
	rm -f qr.o qr_gs.o qr_rot.o qr_ref.o qr
