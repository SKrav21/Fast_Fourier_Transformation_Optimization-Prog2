CFLAGS = -O0 -DN=3 -DDEBUG -lm

default: prog2

clean:
	rm -f *.o
	rm -f prog2

prog2: prog2.c
	gcc ${CFLAGS} prog2.c
