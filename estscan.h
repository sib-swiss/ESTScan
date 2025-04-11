/* $Id: estscan.h,v 1.5.2.3 2002/02/20 12:24:14 clottaz Exp $
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
extern double ComputeGC(char *);
extern int StoreTransits(int, int, int, int, int, int, int, int);
extern int Compute(char *, int, int, int, SV *, int, int, int, int, int);
