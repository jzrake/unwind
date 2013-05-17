# ------------------------------------------------------------------------------
#
# This is the unwind code! It might still be called da-code, we're not sure yet
#
# Authors:
#
# Guido, Marge, Roberto, Jonathan
#
# New York University
#
# ------------------------------------------------------------------------------

CFLAGS = -Wall -lm -O3

default : unwind

unwind : unwind.c
	$(CC) $(CFLAGS) -Wall $< -o $@

clean :
	rm -f unwind
