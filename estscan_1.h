/* $Id: estscan_1.h,v 1.1.1.1 2006/12/18 15:44:47 c4chris Exp $
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

extern int CreateMatrix(int, int, int, SV *, int, int);
extern int Compute(char *, int, int, int, int, SV *, int, int, int, int, int, char *);
