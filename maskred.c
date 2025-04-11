/*
 * $Id: maskred.c,v 1.2 2005/03/10 15:08:04 c4chris Exp $
 *
 * maskred.c
 *
 * Reads from stdin in FASTA format expecting nucleotide data, masks
 * reoccuring tuples with 'N' characters if they overlap by a
 * specified number of nucleotides and writes in FASTA format on
 * stdout.
 *
 *    Usage: maskred [-s <tuplesize>] [-o <overlap>] 
 *                   [-m <min-mask>} < infile > outfile
 *
 * '-s' specifies the tuples used to determine reoccurence. Tuples
 * observed in different frames are not considered reoccuring. Tuples
 * observed in UTRs are also considered from those in CDS and vice
 * versa. '-o' specifies how many nucleotides of subsequnt
 * reoccuring tuples must overlap in order to consider them as one
 * region to be masked. Only regions longer than the value specified
 * with '-m' are actually masked.
 *
 * maskred expects CDS annotation in the FASTA-header of the
 * inputs. Immediately after the tag 'CDS: ' the next two'integers
 * separated by a <space> are considered first and last position of
 * the CDS. CDS is not recognized correctly if its specification does
 * not entirely occur in the first 1023 bites of the header.
 *
 * written by Claudio Lottaz (SIB-ISREC) in September/October 2001 
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#ifndef __GNUC__
#include <alloca.h>
#endif
#include <math.h>
#include <limits.h>

/*******************************************************************************
 *
 *   Command line arguments 
 */

static void usage(char *arg0)
{
  fprintf(stderr, "Usage: %s [options] < infile > outfile\n"
	  "  where options are:\n"
	  "    -m <int>     minimum length of region to be filtered\n"
	  "    -o <int>     overlap required between following recurring tuples\n"
	  "    -s <int>     size of redundancy filter\n"
	  "    -d           debug\n",
	  arg0);
  exit(1);
}

typedef struct _options_t {
  int tupsize;
  int overlap;
  int minmask;
  int debug;
} options_t;

static options_t options;

static void getOptions(int argc, char *argv[])
{
  /* default values */
  options.tupsize = 12;
  options.overlap =  0;
  options.minmask = 30;
  options.debug = 0;

  /* read command line */
  while (1) {
    int c = getopt(argc, argv, "dhm:o:s:");
    if (c == -1) break;
    switch (c) {
    case 'd': options.debug   = 1;            break;
    case 'm': options.minmask = atoi(optarg); break;
    case 'o': options.overlap = atoi(optarg); break;
    case 's': options.tupsize = atoi(optarg); break;
    case 'h': usage(argv[0]);
    default: printf ("Option -%c is unknown\n", c);
      usage(argv[0]);
    }
  }

  /* check switches */
  if (options.tupsize > 16) {
    fprintf(stderr, "maskred: redundancy filter more than 16 nucleotides wide (%d)\n",
	    options.tupsize);
    exit(1);
  }
  if (options.overlap == 0) options.overlap = options.tupsize - 1;
  if (options.overlap >= options.tupsize) {
    fprintf(stderr, "maskred: too much overlap required (%d >= %d)\n",
	    options.overlap, options.tupsize);
    exit(1);
  }
  if (optind != argc) usage(argv[0]);
}

/*******************************************************************************
 *
 *   various data types and utilities
 */

/* set and test single bits an a large array if bits **************************/

/* Store which tuples have been encountered.
   One array per frame and one for UTRs  */
static unsigned char *seenTuples[4];
static unsigned long tupleIndexMask;
static unsigned char bitValue[] = {1, 2, 4, 8, 16, 32, 64, 128};

inline void setBit(unsigned long index, unsigned char *bitField) 
{ bitField[index >> 3] |= bitValue[index & 7]; }

inline int testBit(unsigned long index, unsigned char *bitField) 
{ return(bitField[index >> 3] & bitValue[index & 7]); }

/* Interpret nucleotides ******************************************************/

static char decode[] = {'A', 'C', 'G', 'T', 'N'};
const unsigned int N = 4;

static unsigned long getCode(int c)
{
  unsigned long i = c & 0x1f; /* Get lower bits, get rid of upper/lower info.  */
  switch (i) {
  case  1: return 0; /* This is A.  */
  case  3: return 1; /* This is C.  */
  case  7: return 2; /* This is G.  */
  case 20: return 3; /* This is T.  */
  case  2 :          /* This is B.  */
  case  4:           /* This is D.  */
  case  8:           /* This is H.  */
  case 11:           /* This is K.  */
  case 13:           /* This is M.  */
  case 14:           /* This is N.  */
  case 18:           /* This is R.  */
  case 19:           /* This is S.  */
  case 22:           /* This is V.  */
  case 23:           /* This is W.  */
  case 25: return N; /* This is Y.  */
  default: /* Everything else.  */
    fprintf(stderr, "Bad character in getCode: %c(%d)\n", (char)c, c);
    exit(1);
  }
}

/********************************************************************************
 *
 *   Data analysis
 */

static void printTuple(unsigned long index, int len)
{
  int i;
  unsigned char *t = (unsigned char *)alloca(sizeof(unsigned char) * len);
  for (i = len-1; i >= 0; i--) { t[i] = index & 3; index >>= 2; }
  for (i = 0; i < len; i++) putchar(decode[t[i]]);     
}

inline char get_nucleotide() 
{
  char c = getchar();
  while ((c == 12) || (c == 10)) c = getchar(); /* skip <CR> and <LF> */
  return c;
}

int maskrun(char *string, int *begin, int *end)
     /* in 'string' fills substring from 'begin' to 'end' with Ns.
      * perform this action only of the run is long enough.
      */
{
  int i, len = 0;
  if ((*end - *begin) > options.minmask) {
    for(i = *begin; i <= *end; i++) string[i] = 'N'; 
    len = *end - *begin + 1;
  }
  *begin = *end = -1;
  return len;
}

void mask_file()
{
  int i, pos, size, skip, masked = 0;
  int tupleIndex, maskStart, maskEnd;
  int cdsStart, cdsEnd;
  char c, buf[1024], *longbuf, *s;

  size = 1024;
  longbuf = malloc(sizeof(char)*size);

  c = getchar();
  while(c == '>') {
    /* read and print the header, find CDS */
    fgets(buf, 1024, stdin); 
    fprintf(stdout, ">%s", buf);
    s = strstr(buf, "CDS: ") + 5;
    cdsStart = atoi(s);
    s = strchr(s, ' ') + 1;
    cdsEnd = atoi(s);
    while(buf[strlen(buf)-1] != '\n') { 
      fgets(buf, 1024, stdin); 
      fprintf(stdout, "%s", buf);
    }
    
    /* read sequence into the buffer */
    pos = 0; 
    c = get_nucleotide();
    while((c != '>') && (c != EOF)) {
      if (pos == size) {
	size += 1024;
	longbuf = realloc(longbuf, sizeof(char) * size); 
      }
      longbuf[pos] = c; pos++;
      c = get_nucleotide();
    }

    /* mask redundant tuples */
    tupleIndex = 0;
    maskStart = maskEnd = -1; /* -1 means current region unmasked */
    skip = options.tupsize;   /* tupleIndex is invalid for the next skip positions */
    for (i = 0; i < pos; i++) {
      unsigned long code = getCode(longbuf[i]);
      if (code < N) tupleIndex = ((tupleIndex << 2) | code) & tupleIndexMask; 
      else {
	skip = options.tupsize; 
	maskrun(longbuf, &maskStart, &maskEnd);
      }
      if (skip) skip--;
      else {
	int f = ((cdsStart <= i) && (i <= cdsEnd)) ? (i - cdsStart) % 3 : 3;
	if (testBit(tupleIndex, seenTuples[f])) {
	  if (maskStart == -1) maskStart = i - options.tupsize + 1; 
	  maskEnd = i;
	}
	else {
	  if ((maskStart != -1) && (i - maskEnd) > (options.tupsize - options.overlap)) 
	    masked += maskrun(longbuf, &maskStart, &maskEnd);
	  setBit(tupleIndex, seenTuples[f]);
	}
      }
    }
    masked += maskrun(longbuf, &maskStart, &maskEnd);    

    /* write masked buffer to output */
    for (i = 0; i < pos; i++) {
      putchar(longbuf[i]);
      if (((i+1) % 80) == 0) putchar('\n');
    }
    putchar('\n');
  }
  fprintf(stdout, ">masked nucleotides: %d\n", masked);
  free(longbuf);
}


/********************************************************************************
 *
 *   Main
 */

int main(int argc, char *argv[])
{
  int i;
  unsigned long storeSize; 

  getOptions(argc, argv);

  /* reserve one bit for each nucleotide tuple of the size of the filter */
  storeSize = 1<<(2*options.tupsize - 3);
  for (i = 0; i < 4; i++) {
    seenTuples[i] = (unsigned char *)malloc(sizeof(unsigned char) * storeSize);
    memset(seenTuples[i], 0, storeSize);
  }
  tupleIndexMask = ((storeSize - 1) << 3) | 7 ;
  mask_file();

  /* clean up */
  for (i = 0; i < 4; i++) { free(seenTuples[i]); }
  return 0;
}

/*
 * End of File
 *
 *******************************************************************************/
