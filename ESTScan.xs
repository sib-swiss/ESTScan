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
CreateMatrix(a, b, c, d, e, f)
	int a
	int b
	int c
	SV *d
	int e
	int f
	OUTPUT:
	RETVAL

int
Compute(a, b, c, d, e, f, g, h, i, j, k, l)
	char *a
	int b
	int c
	int d
	int e
	SV *f
	int g
	int h
	int i
        int j
        int k
 	char *l
	OUTPUT:
	RETVAL

# $Id: ESTScan.xs,v 1.7 2000/09/06 07:43:19 clottaz Exp $
#
# Christian Iseli, LICR ITO, Christian.Iseli@licr.org
#
# Copyright (c) 1999 Swiss Institute of Bioinformatics. All rights reserved.











