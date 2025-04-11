# $Id: Makefile.PL,v 1.1.2.2 2002/02/25 13:39:37 clottaz Exp $
#
# Claudio Lottaz, SIB-ISREC
#
# Copyright (c) 2002 Swiss Institute of Bioinformatics. All rights reserved.

print "Writing Makefile for MkTables\n";
open(FH, ">Makefile");
print FH <<EndOfBuffer;
CFLAGS = -O2
CLIBS = -lm

all:	makesmat maskred

maskred: maskred.c
	gcc \$(CFLAGS) -o maskred maskred.c \$(CLIBS)

makesmat: makesmat.c
	gcc \$(CFLAGS) -o makesmat makesmat.c \$(CLIBS)

clean:
	rm Makefile maskred makesmat

realclean:
	rm Makefile maskred makesmat

test : ;
install : ;

EndOfBuffer
close FH;
