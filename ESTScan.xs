#ifdef __cplusplus
extern "C" {
#endif
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#ifdef __cplusplus
}
#endif

#include "estscan.h"


MODULE = ESTScan		PACKAGE = ESTScan
PROTOTYPES: ENABLE

# second parameter is a Perl Array Ref
int
CreateMatrix(a, b, c, d, e, f, g, h, i)
	int a
	int b
	int c
	int d
	int e
	int f
	SV *g
	int h
	int i
	OUTPUT:
	RETVAL

int
StoreTransits(a, b, c, d, e, f, g, h)
	int a
	int b
	int c
	int d
	int e
	int f
	int g
	int h
	OUTPUT:
	RETVAL

double
ComputeGC(a)
	char *a
	OUTPUT:
	RETVAL

int
Compute(a, b, c, d, e, f, g, h, i, j)
	char *a
	int b
	int c
	int d
	SV *e
	int f
	int g
	int h
        int i
	int j
	OUTPUT:
	RETVAL

# $Id: ESTScan.xs,v 1.7.2.3 2002/02/20 12:24:11 clottaz Exp $
#
# Christian Iseli, LICR ITO, Christian.Iseli@licr.org
#
# Copyright (c) 1999 Swiss Institute of Bioinformatics. All rights reserved.
