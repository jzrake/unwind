
CFLAGS = -Wall -lm -O3

default : unwind

unwind : unwind.c
	$(CC) $(CFLAGS) -Wall $< -o $@

clean :
	rm -f unwind

