all:	
	make fin
	make rod
#
# fin
# ---
#
fin: fin.o
	f77 -o fin fin.o
#
# rod
# ---
#
rod: rod.o
	f77 -o rod rod.o
#
# clean
# ---
#
clean:
	rm -f fin.o fin
	rm -f rod.o rod rod.out
