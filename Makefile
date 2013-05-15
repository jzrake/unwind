
CFLAGS = -Wall -lm -O3

default : ma-code

ma-code : ma-code.c
	$(CC) $(CFLAGS) -Wall $< -o $@

clean :
	rm -f ma-code

