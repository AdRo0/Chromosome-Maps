PROG = main

CC = gcc

CFLAGS = `pkg-config --cflags gtk+-3.0`

LIBS = `pkg-config --libs gtk+-3.0` -lm

${PROG}: ${PROG}.c
	${CC} -o ${PROG} ${PROG}.c -Wall ${CFLAGS} ${LIBS} -export-dynamic

clean:
	rm -f ${PROG}
