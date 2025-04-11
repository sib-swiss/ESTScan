/* $Id: estscan_.h,v 1.2 2006/12/14 23:59:58 c4chris Exp $
 *
 * Christian Iseli, LICR ITO, Christian.Iseli@licr.org
 *
 * Copyright (c) 1999 Swiss Institute of Bioinformatics. All rights reserved.
 */
/* Interface file for our EST scanning library stuff.  */

#include <stdlib.h>
#include <time.h>
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

extern int CreateMatrix(int, int, int, int, int, int, SV *, int, int);
extern double ComputeGC(const char *);
extern int StoreTransits(int, int, int, int, int, int, int, int);
extern int Compute(const char *, int, int, int, SV *,
		   int, int, int, int, int);
