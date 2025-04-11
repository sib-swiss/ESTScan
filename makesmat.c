/*
 * $Id: makesmat.c,v 1.4 2005/04/06 23:42:52 c4chris Exp $
 *
 * makesmat.c
 *
 * Reads from stdin in FASTA format expecting full-length messenger
 * RNA data, counts tuples which do not contain ambiguous codes and
 * deduces log-odd emission probabilities for untranslated region and
 * coding sequence as well as positionspecific scoring matrices for
 * start and stop sites. Output is generated in GENSCAN format.
 *
 *    Usage: maskred [-t <tuplesize>] [-p pseudocounts] 
 *                   [-o <s|p|c>] < infile > outfile
 *
 * '-t' specifies the tuples used to determine counts. Tuples observed
 * in different frames are counted apart. Tuples observed in UTRs are
 * also considered separate from those in CDS. '-p' specifies the
 * pseudocounts added in the end. Pseudocounts are added proportional
 * to multiplied single nucleotide occurence and sum up to the number
 * specified using '-p'. The option '-o' allows to select whether
 * counts (c), probabilities (p) or log-odd scores (s) are computed.
 *
 * makersmat expects CDS annotation in the FASTA-header of the
 * inputs. Immediately after the tag 'CDS: ' the next two'integers
 * separated by a <space> are interpreted as the first and last
 * position of the CDS. CDS is not recognized correctly if its
 * specification does not entirely occur in the first 1023 bites of
 * the header.
 *
 * written by Claudio Lottaz (SIB-ISREC) in October 2001 
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

typedef struct _options_t {
  double scoreFactor;
  int tupsize;
  int pseudocounts;
  int startFrames;
  int startOffset;
  int stopFrames;
  int stopOffset;
  char output_type; /* 's':scores (default), 'p':probabilities, 'c':counts */
  int minscore;
  int debug;
} options_t;

static options_t options;

/*******************************************************************************
 *
 *   Command line arguments 
 */

static void
usage(char *arg0)
{
  fprintf(stderr, "Usage: %s [options] < infile > outfile\n"
	  "  where options are:\n"
	  "    -t <int>   tuple size [%d]\n"
	  "    -p <int>   pseudocounts to be added [%d]\n"
	  "    -f <int>   number of frames in start profiles [%d]\n"
	  "    -o <int>   site offset within start profiles [%d]\n"
	  "    -F <int>   number of frames in stop profiles [%d]\n"
	  "    -O <int>   site offset within stop profiles [%d]\n"
	  "    -T <s|p|c> output type: scores, probabilities or counts [%c]\n"
	  "    -m <int>   minimum score [%d]\n"
	  "    -s <float> score multiplication factor [%.1f]\n"
	  "    -h         display usage info\n"
	  "    -d         debug\n",
	  arg0, options.tupsize, options.pseudocounts, options.startFrames,
	  options.startOffset, options.stopFrames, options.stopOffset,
	  options.output_type, options.minscore, options.scoreFactor);
  exit(1);
}

static void
getOptions(int argc, char *argv[])
{
  /* default values */
  options.scoreFactor = 5.0;
  options.tupsize = 6;
  options.pseudocounts =  1;
  options.output_type = 's';
  options.startFrames = 18;
  options.startOffset = 7;
  options.stopFrames = 18;
  options.stopOffset = 6;
  options.minscore = -100;
  options.debug = 0;

  /* read command line */
  while (1) {
    int c = getopt(argc, argv, "t:p:f:F:o:O:T:m:s:dh");
    if (c == -1) break;
    switch (c) {
    case 't': options.tupsize      = atoi(optarg); break;
    case 'p': options.pseudocounts = atoi(optarg); break;
    case 'f': options.startFrames  = atoi(optarg); break;
    case 'o': options.startOffset  = atoi(optarg); break;
    case 'F': options.stopFrames   = atoi(optarg); break;
    case 'O': options.stopOffset   = atoi(optarg); break;
    case 'T': options.output_type  = optarg[0];    break;
    case 'm': options.minscore     = atoi(optarg); break;
    case 's': options.scoreFactor  = atof(optarg); break;
    case 'd': options.debug   = 1;                 break;
    case 'h': usage(argv[0]);
    default: printf ("Option -%c is unknown\n", c);
      usage(argv[0]);
    }
  }

  /* check switches */
  if (options.tupsize > 16) {
    fprintf(stderr, "makesmat: tuplesize too large (%d > 16)\n",
	    options.tupsize);
    exit(1);
  }
  if (options.tupsize < 2) {
    fprintf(stderr, "makesmat: tuplesize too small (%d < 2)\n",
	    options.tupsize);
    exit(1);
  }
  if (options.pseudocounts < 0) {
    fprintf(stderr, "makesmat: negative pseudocounts (%d)\n",
	    options.pseudocounts);
    exit(1);
  }
  if ((options.startFrames - options.startOffset + 1)< options.tupsize) {
    fprintf(stderr,
	    "makesmat: start profile, tuple size to large (%d-%d>=%d)\n", 
	    options.startFrames, options.startOffset, options.tupsize);
    exit(1);
  }
  if ((options.stopFrames - options.stopOffset) < options.tupsize) {
    fprintf(stderr, "makesmat: stop profile, tuple size to large (%d-%d>%d)\n", 
	    options.stopFrames, options.stopOffset, options.tupsize);
    exit(1);
  }
  if ((options.output_type != 's')
      && (options.output_type != 'p')
      && (options.output_type != 'c')) {
    fprintf(stderr, "makesmat: unrecognized output-type (%c)\n",
	    options.output_type);
    exit(1);
  }
  if (optind != argc)
    usage(argv[0]);
}

/******************************************************************************
 *
 *   various data types and utilities
 */

/* Interpret nucleotides  */

static char decode[] = {'A', 'C', 'G', 'T', 'N'};
const int N = 4;

static char
getCode(int c)
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

static inline char
get_nucleotide(void) 
{
  char c = getchar();
  while ((c == 12) || (c == 10)) c = getchar(); /* skip <CR> and <LF> */
  return c;
}

static void
print_tuple(unsigned long index, int len)
{
  int i;
  unsigned char *t = (unsigned char *)alloca(sizeof(unsigned char) * len);
  for (i = len-1; i >= 0; i--) {
    t[i] = index & 3;
    index >>= 2;
  }
  for (i = 0; i < len; i++)
    putchar(decode[t[i]]);     
}

/* tuple and nucleotide counters  */

static unsigned long singleTotal;  /* number of nucleotides */
static unsigned long singleCtr[4]; /* occurance of single nucleotides */
static unsigned long *startctr[4]; /* counters for start PSSM */
static unsigned long *stopctr[4];  /* counters for stop PSSM */

typedef struct _counters_t {
  unsigned long tupsize1Total;     /* number of (tupsize-1)-tuples */
  unsigned long tupsizeTotal;      /* number of (tupsize)-tuples */
  unsigned long *tupsize1Ctr;      /* occurence of (tupsize-1)-tuples */
  unsigned long *tupsizeCtr;       /* occurence of tupsize-tuples */
} counters_t, *counters_p_t;
static counters_t ctr[4];          /* counters for frames 0, 1 and 2 as well as UTR */

static void
initCounters(void) 
{
  int i;
  int nbTuples1 = (1 << (2*(options.tupsize - 1)));
  int nbTuples =  (1 << (2*options.tupsize));
  for (i = 0; i < 4; i++) {
    ctr[i].tupsize1Total = 0;
    ctr[i].tupsizeTotal = 0;
    ctr[i].tupsize1Ctr
      = (unsigned long *)malloc(sizeof(unsigned long) * nbTuples1);
    ctr[i].tupsizeCtr
      = (unsigned long *)malloc(sizeof(unsigned long) * nbTuples);
    memset(ctr[i].tupsize1Ctr, 0, sizeof(unsigned long) * nbTuples1);
    memset(ctr[i].tupsizeCtr,  0, sizeof(unsigned long) * nbTuples);

    startctr[i]
      = (unsigned long *)malloc(sizeof(unsigned long) * options.startFrames);
    memset(startctr[i], 0, sizeof(unsigned long) * options.startFrames);
    stopctr[i]
      = (unsigned long *)malloc(sizeof(unsigned long) * options.stopFrames);
    memset(stopctr[i], 0, sizeof(unsigned long) * options.stopFrames);
  }

  singleTotal = 0;
  memset(singleCtr, 0, sizeof(unsigned long) * 4);
}

static void
update_counters(unsigned long index, int frame, int skip) 
/* skip indicates how many high order nts in index are not valid */
{
  if (skip <= 1) { 
    ctr[frame].tupsize1Ctr[index >> 2]++; 
    ctr[frame].tupsize1Total++;
  }
  if (skip == 0) { 
      ctr[frame].tupsizeCtr[index]++; 
      ctr[frame].tupsizeTotal++;
    }
}

static void
free_counters(void)
{
  int i;
  for (i = 0; i <= 3; i++) { 
    free(ctr[i].tupsize1Ctr); 
    free(ctr[i].tupsizeCtr); 
  }
}

 /******************************************************************************
  *
  *   Data analysis
  */

void
count_tuples(void)
{
  int i, j, pos, size, skip;
  int tupleIndex,  tupleIndexMask;
  int cdsStart, cdsEnd;
  char c, buf[1024], *s;
  unsigned char *longbuf;

  size = 1024;
  longbuf = malloc(sizeof(unsigned char)*size);
  tupleIndexMask = (1<<(2*options.tupsize)) - 1;

  c = getchar();
  while (c == '>') {
    /* read and print the header, find CDS */
    fgets(buf, 1024, stdin); 
    s = strstr(buf, "CDS: ");
    if (s == NULL) { /* skip to next entry */
      c = getchar();
      while((c != EOF) && (c != '>'))
	c = getchar();
      continue;
    }
    s += 5;
    cdsStart = atoi(s) - 1; /* C index start at 0! */
    s = strchr(s, ' ') + 1;
    cdsEnd =   atoi(s) - 1; /* C index start at 0! */
    while (buf[strlen(buf)-1] != '\n') {
      fgets(buf, 1024, stdin);
    }

    /* read sequence into the buffer */
    pos = 0; 
    c = get_nucleotide();
    while ((c != '>') && (c != EOF)) {
      if (pos == size) {
	size += 1024;
	longbuf = realloc(longbuf, sizeof(unsigned char) * size); 
      }
      longbuf[pos] = getCode(c);
      pos++;
      c = get_nucleotide();
    }

     /* count single nucleotides */
    for (i = 0; i < pos; i++) {
      if (longbuf[i] < N) {
	singleCtr[longbuf[i]]++;
	singleTotal++;
      }
    }

    /* count 5'UTR */
    tupleIndex = 0;
    /* tupleIndex is invalid for the next skip positions */
    skip = options.tupsize;
    for (i = 0; i <= cdsStart - options.startOffset; i++) {
      if (longbuf[i] < N)
	tupleIndex = ((tupleIndex << 2) | longbuf[i]) & tupleIndexMask;
      else
	skip = options.tupsize; 
      if (skip)
	skip--;
      update_counters(tupleIndex, 3, skip); /* update UTR counters */
    }

    /* count start profile */
    j = cdsStart - options.startOffset + 1;
    i = (j < 0) ? 0 : j;
    while (i < j + options.startFrames) {
      if (longbuf[i] < N)
	startctr[longbuf[i]][i - j]++; 
      i++;
    }

    /* count CDS core */
    tupleIndex = 0;
    /* tupleIndex is invalid for the next skip positions */
    skip = options.tupsize;
    for (i = cdsStart + options.startFrames - options.startOffset + 1; 
	 i <= cdsEnd - options.stopOffset;
	 i++) {
      if (longbuf[i] < N)
	tupleIndex = ((tupleIndex<<2) | longbuf[i]) & tupleIndexMask;
      else
	skip = options.tupsize; 
      if (skip)
	skip--;
      /* update UTR counters */
      update_counters(tupleIndex, (i - cdsStart)%3, skip);
    }

    /* count stop profile */
    for (i = cdsEnd - options.stopOffset + 1; 
	 (i <= cdsEnd + options.stopFrames - options.stopOffset) && (i < pos);
	 i++) {
      if (longbuf[i] < N)
	stopctr[longbuf[i]][i - cdsEnd + options.stopOffset - 1]++;
    }

    /* count 3'UTR */
    tupleIndex = 0;
    /* tupleIndex is invalid for the next skip positions */
    skip = options.tupsize;
    for (i = cdsEnd + options.stopFrames - options.stopOffset + 1;
	 i < pos;
	 i++) {
      if (longbuf[i] < N)
	tupleIndex = ((tupleIndex<<2) | longbuf[i]) & tupleIndexMask; 
      else
	skip = options.tupsize; 
      if (skip)
	skip--;
      update_counters(tupleIndex, 3, skip); /* update UTR counters */
    }
  }
  free(longbuf);
}

 /****************************************************************************
  *
  *   Print tables
  */

static double
pseudo(int tuple, int tupsize)
{
  int i;
  double p = 1.0;
  for (i = 0; i < tupsize; i++) {
    p *= ((double) singleCtr[tuple & 3]) / ((double) singleTotal);
    tuple >>= 2;
  }
  return p * options.pseudocounts;
}

static void
print_cdstable(double *singleProb) 
{
  int f, i, *scores[3];
  double *tuple1Prob[3], *tupleProb[3];

  int nbTuples1 = (1 << 2*(options.tupsize-1));
  int nbTuples = (1 << 2*options.tupsize);

  /* compute probabilities */
  for (f = 0; f < 3; f++) {
    if (options.debug)
      printf("Probabilities for frame %d %d\n", f, options.tupsize - 1);
    tuple1Prob[f] = (double *) malloc(sizeof(double) * nbTuples1);
    for (i = 0; i < nbTuples1; i++) {
      tuple1Prob[f][i] = ((double) ctr[f].tupsize1Ctr[i]
			  + pseudo(i, options.tupsize - 1))
			 / (double) (ctr[f].tupsize1Total
				     + options.pseudocounts); 
      if (options.debug) {
	print_tuple(i, options.tupsize-1);
	printf(":%8g=(%8ld + %8g)/(%8ld + %8d)\n",
	       tuple1Prob[f][i], ctr[f].tupsize1Ctr[i],
	       pseudo(i, options.tupsize-1),
	       ctr[f].tupsize1Total,
	       options.pseudocounts);
      }
    }
    if (options.debug)
      printf("Probabilities for frame %d %d\n", f, options.tupsize);
    tupleProb[f] = (double *)malloc(sizeof(double) * nbTuples);
    for (i = 0; i < nbTuples; i++) {
      tupleProb[f][i] = ((double) ctr[f].tupsizeCtr[i] +
			 pseudo(i, options.tupsize))
			/ (double) (ctr[f].tupsizeTotal
				    + options.pseudocounts); 
      if (options.debug) {
	print_tuple(i, options.tupsize);
	printf(":%8g=(%8ld + %8g) / (%8ld + %8d)\n",
	       tupleProb[f][i],
	       ctr[f].tupsizeCtr[i],
	       pseudo(i, options.tupsize),
	       ctr[f].tupsizeTotal,
	       options.pseudocounts);
      }
    }
  }

  /* compute log-odds and scores */
  for (f = 0; f < 3; f++) { 
    if (options.debug)
      printf("Log-odds for frame %d\n", f);
    scores[f] = (int *)malloc(sizeof(int) * nbTuples);
    for (i = 0; i < nbTuples; i++) {
      if ((tupleProb[f][i] == 0.0) || (tuple1Prob[f][i >> 2] == 0.0)) 
	scores[f][i] = options.minscore;
      else {
	double score = log(tupleProb[f][i]
			   / singleProb[i & 3]
			   / tuple1Prob[f][i >> 2])
		       / M_LN2 * options.scoreFactor;
	scores[f][i]
	  = (score < options.minscore) ? options.minscore : round(score);
      }
      if (options.debug) {
	print_tuple(i, options.tupsize);
	printf(": %8g = %8g / %8g, %4d = 10*log(%8g / %8g)\n", 
	       tupleProb[f][i] / tuple1Prob[f][i >> 2], 
	       tupleProb[f][i],
	       tuple1Prob[f][i >> 2],
	       scores[f][i],
	       tupleProb[f][i] / tuple1Prob[f][i >> 2],
	       singleProb[i & 3]);
      }
    }
  }

  /* print table */
  if (options.debug) { 
    if (options.output_type == 'p') 
      printf("single: %.6f %.6f %.6f %.6f\n",
	     singleProb[0],
	     singleProb[1], 
	     singleProb[2],
	     singleProb[3]);
    if (options.output_type == 'c') 
      printf("single (total %ld): %-6ld %-6ld %-6ld %-6ld\n",
	     singleTotal,
	     singleCtr[0],
	     singleCtr[1],
	     singleCtr[2],
	     singleCtr[3]);
  }
  for (f = 0; f < 3; f++) {
    for (i = 0; i < 1 << (2*options.tupsize); i += 4) {
      int rowIndex = i >> 2;
      if (options.debug) { 
	print_tuple(rowIndex, options.tupsize-1); 
	if (options.output_type == 'p')
	  printf(" (total %.5f)", tuple1Prob[f][i >> 2]); 
	if (options.output_type == 'c')
	  printf(" (total %ld)",ctr[f].tupsize1Ctr[i>>2]);
	printf(": ");
      }
      switch (options.output_type) {
      case 's': 
	printf("%-6d %-6d %-6d %-6d\n",
	       scores[f][i],
	       scores[f][i + 1], 
	       scores[f][i + 2],
	       scores[f][i + 3]);
	break;
      case 'p': 
	printf("%.6f %.6f %.6f %.6f\n",
	       tupleProb[f][i],
	       tupleProb[f][i + 1],
	       tupleProb[f][i + 2],
	       tupleProb[f][i + 3]);
	break;
      case 'c': 
	printf("%-6ld %-6ld %-6ld %-6ld\n",
	       ctr[f].tupsizeCtr[i],
	       ctr[f].tupsizeCtr[i + 1], 
	       ctr[f].tupsizeCtr[i + 2],
	       ctr[f].tupsizeCtr[i + 3]);
      }
    }
  }

  /* clean up */
  for (f = 0; f < 3; f++) { 
    free(tuple1Prob[f]);
    free(tupleProb[f]);
    free(scores[f]);
  }
}

static void
print_utrtable(double *singleProb)
{
  int i, *scores;
  double *tuple1Prob, *tupleProb, currTotal;

  int nbTuples1 = (1 << 2*(options.tupsize-1));
  int nbTuples = (1 << 2*options.tupsize);

  /* compute probabilities */
  currTotal = ctr[3].tupsize1Total;
  tuple1Prob = (double *) malloc(sizeof(double) * nbTuples1);
  for (i = 0; i < nbTuples1; i++) 
    tuple1Prob[i] = ((double) ctr[3].tupsize1Ctr[i]
		     + pseudo(i, options.tupsize-1))
		    / (currTotal + options.pseudocounts); 
  currTotal = ctr[3].tupsizeTotal;
  tupleProb = (double *) malloc(sizeof(double) * nbTuples);
  for (i = 0; i < nbTuples; i++) 
    tupleProb[i] = ((double) ctr[3].tupsizeCtr[i]
		    + pseudo(i, options.tupsize))
		   / (currTotal + options.pseudocounts); 

  /* compute log-odds and scores */
  scores = (int *) malloc(sizeof(int) * nbTuples);
  for (i = 0; i < nbTuples; i++) {
    if ((tupleProb[i] == 0.0) || (tuple1Prob[i>>2] == 0.0))
      scores[i] = options.minscore;
    else {
      double score = log(tupleProb[i] / singleProb[i&3] / tuple1Prob[i>>2])
		     / M_LN2 * options.scoreFactor;
      scores[i] = (score < options.minscore) ? options.minscore : round(score);
    }
  }

  /* print table */
  if (options.debug) { 
    if (options.output_type == 'p') 
      printf("single: %.6f %.6f %.6f %.6f\n",
	     singleProb[0],
	     singleProb[1], 
	     singleProb[2],
	     singleProb[3]);
    if (options.output_type == 'c') 
      printf("single (total %ld): %-6ld %-6ld %-6ld %-6ld\n",
	     singleTotal,
	     singleCtr[0],
	     singleCtr[1],
	     singleCtr[2],
	     singleCtr[3]);
  }
  for (i = 0; i < 1 << (2*options.tupsize); i += 4) {
    int rowIndex = i >> 2;
    if (options.debug) { 
      print_tuple(rowIndex, options.tupsize-1); 
      if (options.output_type == 'p')
	printf(" (total %.5f)", tuple1Prob[i >> 2]); 
      if (options.output_type == 'c')
	printf(" (total %ld)", ctr[3].tupsize1Ctr[i>>2]); 
      printf(": ");
    }
    switch (options.output_type) {
    case 's': 
      printf("%-6d %-6d %-6d %-6d\n",
	     scores[i],
	     scores[i + 1], 
	     scores[i + 2],
	     scores[i + 3]);
      break;
    case 'p': 
      printf("%.6f %.6f %.6f %.6f\n",
	     tupleProb[i],
	     tupleProb[i + 1], 
	     tupleProb[i + 2],
	     tupleProb[i + 3]);
      break;
    case 'c': 
      printf("%-6ld %-6ld %-6ld %-6ld\n",
	     ctr[3].tupsizeCtr[i],
	     ctr[3].tupsizeCtr[i + 1], 
	     ctr[3].tupsizeCtr[i + 2],
	     ctr[3].tupsizeCtr[i + 3]);
    }
  }

  /* clean up */
  free(tuple1Prob);
  free(tupleProb);
  free(scores);
}

static void
print_pssm(double *singleProb, unsigned long **ctr, int frames)
{
  int i, j;
  for (i = 0; i < frames; i++) {
    int currTotal = 0;
    for (j = 0; j < 4; j++) {
      currTotal += ctr[j][i];
    }
    for (j = 0; j < 4; j++) {
      double pseudo = options.pseudocounts
		      * (double) singleCtr[j] / (double) singleTotal;
      double p = ((double) ctr[j][i] + pseudo)
		 / ((double) currTotal + options.pseudocounts);
      switch(options.output_type) {
      case 'c':
	printf("%-6ld ", ctr[j][i]);
	break;
      case 'p':
	printf("%.6f ", p);
	break;
      case 's': 
	if ((p == 0.0) || (singleProb[j] == 0.0))
	  printf("%-6d ", options.minscore);
	else {
	  double score = log(p / singleProb[j]) / M_LN2 * 10.0;
	  printf("%-6d ",
		 (score < options.minscore)
		  ? options.minscore
		  : (int) round(score));
	}
      }
    }
    printf("\n");
  }
}

/*****************************************************************************
 *
 *   Main
 */

int
main(int argc, char *argv[])
{
  int i;
  double currTotal, singleProb[4];

  /* initialize and count tuples */
  getOptions(argc, argv);
  initCounters();
  count_tuples();

  /* compute probabilities of single nucleotides */
  currTotal = singleTotal;
  for (i = 0; i < 4; i++)
    singleProb[i] = (double) singleCtr[i] / currTotal; 

  /* output tables */
  printf("FORMAT: <NAME> CODING REGION %d 3 1 %c C+G: <CG>\n",
	 options.tupsize, options.output_type);
  print_cdstable(singleProb);

  printf("FORMAT: <NAME> UNTRANSLATED REGION %d 1 1 %c C+G: <CG>\n",
	 options.tupsize, options.output_type);
  print_utrtable(singleProb);

  printf("FORMAT: <NAME> START PROFILE 1 %d %d %c C+G: <CG>\n",
	 options.startFrames, options.startOffset, options.output_type);
  print_pssm(singleProb, startctr, options.startFrames);

  printf("FORMAT: <NAME> STOP PROFILE 1 %d %d %c C+G: <CG>\n",
	 options.stopFrames, options.stopOffset, options.output_type);
  print_pssm(singleProb, stopctr, options.stopFrames);

  /* clean up */
  free_counters();
  return 0;
}

/*
 * End of File
 *
 ****************************************************************************/
