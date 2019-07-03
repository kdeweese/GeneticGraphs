CC = g++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
CPPFLAGS = -std=c++11 $(DEFS)
LFLAGS = -Wall $(DEBUG)

all : tester create_random mutate recombine mutateandrecombine intersection attach

tester: tester.o
	$(CC) $(LFLAGS) tester.o -o tester

tester.o: tester.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) tester.cpp

attach: attach.o
	$(CC) $(LFLAGS) attach.o -o attach

attach.o: attach.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) attach.cpp

create_random: create_random.o
	$(CC) $(LFLAGS) create_random.o -o create_random

create_random.o: create_random.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) create_random.cpp

mutate: mutate.o
	$(CC) $(LFLAGS) mutate.o -o mutate

mutate.o: mutate.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) mutate.cpp

recombine: recombine.o
	$(CC) $(LFLAGS) recombine.o -o recombine

recombine.o: recombine.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) recombine.cpp

mutateandrecombine: mutateandrecombine.o
	$(CC) $(LFLAGS) mutateandrecombine.o -o mutateandrecombine

mutateandrecombine.o: mutateandrecombine.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) mutateandrecombine.cpp

intersection: intersection.o
	$(CC) $(LFLAGS) intersection.o -o intersection

intersection.o: intersection.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) intersection.cpp


clean:
	\rm *.o tester create_random mutate recombine mutateandrecombine intersection attach
