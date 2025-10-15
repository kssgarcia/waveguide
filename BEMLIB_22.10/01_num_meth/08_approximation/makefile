#
# bernstein
# ---
#
bernstein: bernstein.o 
	f77 -o bernstein bernstein.o 
#
# bspl_un
# ---
#
bspl_un: bspl_un.o 
	f77 -o bspl_un bspl_un.o 
#
# bspline
# -------
#
bspline: bspline.o 
	f77 -o bspline bspline.o
#
# fourier_id
# ---
#
fourier_id: fourier_id.o 
	f77 -o fourier_id fourier_id.o 
#
# fft
# ---
#
fft: fft_dr.o fft.o
	f77 -o fft fft_dr.o fft.o
fft.o: fft.f
	f77 -c fft.f
#
# ft
# --
#
ft: ft.o
	f77 -o ft ft.o 
ft.o: ft.f
	f77 -c ft.f
#
# ft_2d
# -----
#
ft_2d: ft_2d.f
	f77 -o ft_2d ft_2d.f
#
# all
# ---
#
all:
	make bernstein
	make bspline
	make bspl_un
	make fourier_id
	make fft
	make ft
	make ft_2d
#
# clean
# -----
#
clean:
	rm -f core PLOTDAT
	rm -f bernstein.o bernstein
	rm -f bspl_un.o bspl_un
	rm -f bspline.o bspline
	rm -f fourier_id.o fourier_id
	rm -f fft_dr.o fft.o fft fft.out
	rm -f ft.o ft ft.out
	rm -f ft_2d.o ft_2d
#
# purge
# -----
#
purge:
	rm -f core PLOTDAT
	rm -f bernstein.o 
	rm -f bspl_un.o
	rm -f bspline.o
	rm -f fourier_id.o
	rm -f fft_dr.o fft.o
	rm -f ft.o ft.out
	rm -f ft_2d.o
