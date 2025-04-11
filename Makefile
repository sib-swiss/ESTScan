# $Id: Makefile,v 1.1 2006/05/01 09:22:58 c4chris Exp $
# Set the appropriate compilers and options for your system:
# Any system with GNU compilers:
 CC = gcc
 CFLAGS = -O2
 F77 = g77
 FFLAGS = -O2
 LDFLAGS = -lm

# Linux with Intel compilers:
# CC = icc
# CFLAGS = -O3 -ipo -axP
# F77 = ifort
# FFLAGS = -O3 -ipo -axP

PROGS=maskred makesmat estscan winsegshuffle

all: $(PROGS)

clean:
	\rm -f *~ $(PROGS) *.o

maskred: maskred.o
	$(CC) $(LDFLAGS) -o $@ $<

makesmat: makesmat.o
	$(CC) $(LDFLAGS) -o $@ $<

estscan: estscan.o
	$(CC) $(LDFLAGS) -o $@ $<

winsegshuffle: winsegshuffle.o
	$(F77) $(LDFLAGS) -o $@ $<

.c.o:
	$(CC) $(CFLAGS) -c $<

.f.o:
	$(F77) $(FFLAGS) -c $<
