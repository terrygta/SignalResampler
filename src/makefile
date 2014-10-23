CC = g++
CFLAGS = -Wall -g
LDFLAGS = -L/usr/lib -lboost_program_options

example : example.cpp resample.cpp resample.h upfirdn.h
	${CC} example.cpp resample.cpp  ${CFLAGS} ${LDFLAGS} -o example

.PHONY: clean

clean:
	rm -f *.o example