CC = g++
CFLAGS = -std=c++17 -Wall -g
LDFLAGS = -L/usr/lib

example : example.cpp resample.cpp resample.h upfirdn.h
	${CC} example.cpp resample.cpp  ${CFLAGS} ${LDFLAGS} -o example

.PHONY: clean

clean:
	rm -f *.o example
