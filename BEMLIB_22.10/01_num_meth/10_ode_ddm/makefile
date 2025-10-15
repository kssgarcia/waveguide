#
# fdm
# ---
#
fdm: fdm.o thomas.o
	f77 -o fdm fdm.o thomas.o
thomas.o:
	f77 -c thomas.f
#
# fem
# ---
#
fem: fem.o thomas.o
	f77 -o fem fem.o thomas.o
#
# all
# -----
#
all:
	make fdm
	make fem
#
# clean
# -----
#
clean:
	rm -f fdm.o thomas.o fdm
	rm -f fem.o thomas.o fem fem.out
#
# purge
# -----
#
purge:
	rm -f fdm.o thomas.o
	rm -f fem.o thomas.o fem.out
