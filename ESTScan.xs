#ifdef __cplusplus
extern "C" {
#endif
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#ifdef __cplusplus
}
#endif

#include "estscan_.h"


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
	const char *a
	OUTPUT:
	RETVAL

int
Compute(a, b, c, d, e, f, g, h, i, j)
	const char *a
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

# $Id: ESTScan.xs,v 1.2 2006/12/14 23:59:58 c4chris Exp $
#
# Christian Iseli, LICR ITO, Christian.Iseli@licr.org
#
# Copyright (c) 1999 Swiss Institute of Bioinformatics. All rights reserved.
