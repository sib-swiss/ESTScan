#ifdef __cplusplus
extern "C" {
#endif
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#ifdef __cplusplus
}
#endif

#include "estscan_1.h"


MODULE = ESTScan1		PACKAGE = ESTScan1
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

# $Id: ESTScan1.xs,v 1.1.1.1 2006/12/18 15:44:47 c4chris Exp $
#
# Christian Iseli, LICR ITO, Christian.Iseli@licr.org
#
# Copyright (c) 1999 Swiss Institute of Bioinformatics. All rights reserved.











