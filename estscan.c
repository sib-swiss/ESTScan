/* $Id: estscan.c,v 1.5 2007/02/01 15:15:11 c4chris Exp $
 *
 * Christian Iseli, LICR ITO, Christian.Iseli@licr.org
 *
 * Copyright (c) 2004 Swiss Institute of Bioinformatics.  All rights reserved.
 *
 * Compile with -std=gnu99
 */
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <locale.h>
#ifdef DEBUG
#include <mcheck.h>
#endif
#if !defined(__GNUC__) && defined(sun)
#define inline
#endif

#define BUF_SIZE 4096
#define MT_UNKNOWN -1
#define MT_CODING 0
#define MT_UNTRANSLATED 1
#define MT_START 2
#define MT_STOP 3
#define MT_COUNT 4

#define min(x, y)       ((x > y) ? (y) : (x))
#define max(x, y)       ((x < y) ? (y) : (x))

typedef struct _matrix {
  signed char **m;
  char *name;
  char *kind;
  double CGmin;
  double CGmax;
  int matType;
  unsigned int order;
  unsigned int frames;
  int offset;
} matrix_t, *matrix_p_t;

typedef struct _read_buf_t {
  char *line;
  unsigned int lmax;
  unsigned int lc;
  unsigned int ic;
  char in[BUF_SIZE];
} read_buf_t, *read_buf_p_t;

typedef struct _seq_t {
  const char *fName;
  char *header;
  unsigned char *seq;
  double GC_pct;
  read_buf_t rb;
  int fd;
  unsigned int len;
  unsigned int maxHead;
  unsigned int max;
} seq_t, *seq_p_t;

typedef struct _result_t {
  unsigned char *s;
  int score;
  unsigned int start;
  unsigned int stop;
  int reverse;
} result_t, *result_p_t;

typedef union _col_elt_t {
  void **elt;
  matrix_p_t *m;
  result_p_t *r;
} col_elt_t;

typedef struct _col_t {
  col_elt_t e;
  unsigned int size;
  unsigned int nb;
} col_t, *col_p_t;

typedef struct _options_t {
  FILE *out;
  FILE *transl;
  char *matrix;
  double percent;
  double both;
  int min;
  int dPen;
  int iPen;
  int ts5uPen;
  int tscPen;
  int ts3uPen;
  int t5ucPen;
  int t5uePen;
  int tc3uPen;
  int tcePen;
  int t3uePen;
  int Nvalue;
  unsigned int sWidth;
  int all;
  int maxOnly;
  int skipLen;
  int minLen;
  int no_del;
  int single;
} options_t;

static const char Version[] =
"This is ESTScan version 3.0.1.\n"
"Copyright (c) 1999-2007 by the Swiss Institute of Bioinformatics.\n"
"All rights reserved. See the file COPYRIGHT for details.\n";

static const char Usage[] =
"%s [options] [<FASTA file> ...]\n\n"
#ifdef DEBUG
"Debug version\n\n"
#endif
"Available options (default value in braces[]):\n"
"  -a          All in one sequence output\n"
"  -b <float>  only results are shown, which have scores higher than this \n"
"              fraction of the best score [%f].\n"
"  -d <int>    deletion penalty [%d]\n"
"  -h          print this usage information\n"
"  -i <int>    insertion penalty [%d]\n"
"  -l <int>    only results longer than this length are shown [%d]\n"
"  -M <file>   score matrices file ($ESTSCANDIR/Hs.smat)\n"
"              [%s]\n"
"  -m <int>    min value in matrix [%d]\n"
"  -N <int>    how to compute the score of N [%d]\n"
"  -n          remove deleted nucleotides from the output\n"
"  -O          report header information for best match only\n"
"  -o <file>   send output to file.  - means stdout.  If both -t and -o specify\n"
"              stdout, only proteins will be written.\n"
"  -p <float>  GC select correction for score matrices [%f]\n"
"  -S          only analyze positive strand\n"
"  -s <int>    Skip sequences shorter than length [%d]\n"
"  -T <int*>   8 integers used as log-probabilities for transitions,\n"
"              start->5'UTR, start->CDS, start->3'UTR, 5'UTR->CDS,\n"
"              5'UTR->end, CDS->3'UTR, CDS->end, 3'UTR->end\n"
"              [%d, %d, %d, %d, %d, %d, %d, %d]\n"
"  -t <file>   Translate to protein.  - means stdout.\n"
"              will go to the file and the nucleotides will still go to stdout.\n"
"  -v          version information\n"
"  -w <int>    width of the FASTA sequence output [%d]\n";

static options_t options;
static char *argv0;

/* declaration of indexes also used in getFrame */
static int iBegin, i5utr, iStart, iCds, iStop, i3utr;
/* next tsize states implement insertion/deletion after nucleotide in frame index */
static int iInsAfter[3], iDelAfter[3];
/* last tsize states implemented insertion/deletion before nucleotide in frame index */
static int iInsNext[3], iDelNext[3];

static unsigned int maxSize = 0;
static int *V  = NULL;
static int *tr = NULL;

static const unsigned char dna_complement[256] =
  "                                                                "
  " TVGH  CD  M KN   YSA BWXR       tvgh  cd  m kn   ysa bwxr      "
  "                                                                "
  "                                                                ";
/* ................................................................ */
/* @ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~. */
/* ................................................................ */
/* ................................................................ */

#ifdef __GNUC__
static void
fatal(const char *fmt, ...)
     __attribute__ ((format (printf, 1, 2) , __noreturn__));
#endif

static void
fatal(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  fflush(stdout);
  if (argv0) {
    char *p = strrchr(argv0, '/');
    fprintf(stderr, "%s: ", p ? p+1 : argv0);
  }
  vfprintf(stderr, fmt, ap);
  va_end(ap);
#ifdef DEBUG
  abort();
#else
  exit(1);
#endif
}

static int
intCompare(const void *a, const void *b)
{
  int ia = * (int *) a;
  int ib = * (int *) b;
  return ia < ib ? -1 : (ia > ib ? 1 : 0);
}

static void *
xmalloc(size_t size)
{
  void *res = malloc(size);
  if (res == NULL)
    fatal("malloc of %zd failed: %s (%d)\n", size, strerror(errno),
	  errno);
  return res;
}

#if 0
static void *
xcalloc(size_t nmemb, size_t size)
{
  void *res = calloc(nmemb, size);
  if (res == NULL)
    fatal("calloc of %zd, %zd failed: %s (%d)\n", nmemb, size,
	  strerror(errno), errno);
  return res;
}
#endif

static void *
xrealloc(void *ptr, size_t size)
{
  void *res = realloc(ptr, size);
  if (res == NULL)
    fatal("realloc of %p to %zd failed: %s (%d)\n", ptr, size,
	  strerror(errno), errno);
  return res;
}

static void
grow_read_buf(read_buf_p_t b)
{
  b->lmax += BUF_SIZE;
  b->line = xrealloc(b->line, b->lmax * sizeof(char));
}

static char *
shuffle_line(read_buf_p_t b, size_t *cur)
{
  if (b->ic == 0 || *cur >= b->ic)
    return NULL;
  /* Make sure we have enough room in line.  */
  if (b->lmax <= b->lc + (b->ic - *cur))
    grow_read_buf(b);
  while (*cur < b->ic && b->in[*cur] != '\n')
    b->line[b->lc++] = b->in[(*cur)++];
  if (*cur < b->ic) {
    /* Ok, we have our string.  */
    /* Copy the newline.  */
    b->line[b->lc++] = b->in[(*cur)++];
    /* We should be fine, since we read BUF_SIZE -1 at most...  */
    b->line[b->lc] = 0;
    /* Adjust the input buffer.  */
    if (*cur < b->ic) {
      memmove(b->in, b->in + *cur, (b->ic - *cur) * sizeof(char));
      b->ic -= *cur;
    } else
      b->ic = 0;
    *cur = 0;
    return b->line;
  }
  /* Go read some more.  */
  b->ic = 0, *cur = 0;
  return NULL;
}

static char *
read_line_buf(read_buf_p_t b, int fd)
{
  char *s = NULL;
  ssize_t rc;
  size_t cur = 0;
  b->lc = 0;
  if ((s = shuffle_line(b, &cur)) != NULL)
    return s;
  do {
    if ((rc = read(fd, b->in + b->ic, BUF_SIZE - b->ic - 1)) == -1) {
      if (errno != EINTR)
	fatal("Could not read from %d: %s(%d)\n",
	      fd, strerror(errno), errno);
    } else
      b->ic += rc;
    s = shuffle_line(b, &cur);
    if (s == NULL && rc == 0) {
      /* Got to the EOF...  */
      b->line[b->lc] = 0;
      s = b->line;
    }
  } while (s == NULL);
  return s;
}

static void
init_buf(read_buf_p_t b)
{
  b->line = xmalloc(BUF_SIZE * sizeof(char));
  b->lmax = BUF_SIZE;
  b->lc = 0;
  b->ic = 0;
}

static void
free_buf(read_buf_p_t b)
{
  free(b->line);
}

static void
init_seq(const char *fName, seq_p_t sp)
{
  sp->fName = fName;
  sp->header = NULL;
  sp->seq = NULL;
  init_buf(&sp->rb);
  if (fName != NULL) {
    sp->fd = open(fName, O_RDONLY);
    if (sp->fd == -1)
      fatal("Could not open file %s: %s(%d)\n",
	    fName, strerror(errno), errno);
  } else
    sp->fd = 0;
  sp->len = 0;
  sp->maxHead = 0;
  sp->max = 0;
  read_line_buf(&sp->rb, sp->fd);
}

static int
get_next_seq(seq_p_t sp)
{
  const int lenStr = 24;
  unsigned int headerLen;
  char *buf = sp->rb.line;
  int res;
  unsigned int ctr[256], gc, atgc;
  ctr['A'] = ctr['C'] = ctr['G'] = ctr['T'] = 0;
  while (sp->rb.lc > 0 && buf[0] != '>')
    buf = read_line_buf(&sp->rb, sp->fd);
  if (sp->rb.lc == 0)
    return -1;
  /* We have the FASTA header.  */
  if (sp->rb.lc + lenStr + 1 > sp->maxHead) {
    sp->maxHead = sp->rb.lc + lenStr + 1;
    sp->header = (char *) xrealloc(sp->header, sp->maxHead * sizeof(char));
  }
  headerLen = sp->rb.lc;
  memcpy(sp->header, buf, (sp->rb.lc + 1) * sizeof(char));
  sp->len = 0;
  buf = read_line_buf(&sp->rb, sp->fd);
  while (sp->rb.lc > 0 && buf[0] != '>') {
    unsigned char c;
    /* Make sure we have enough room for this additional line.  */
    if (sp->len + sp->rb.lc + 1 > sp->max) {
      sp->max = max(sp->len + sp->rb.lc + 1,
		    sp->max + 0x40000);
      sp->seq = (unsigned char *)
	xrealloc(sp->seq, sp->max * sizeof(unsigned char));
    }
    while ((c = *buf++) != 0) {
      if (isupper(c)) {
	ctr[c] += 1;
	sp->seq[sp->len++] = c;
      } else if (islower(c)) {
	c = toupper(c);
	ctr[c] += 1;
	sp->seq[sp->len++] = c;
      }
    }
    buf = read_line_buf(&sp->rb, sp->fd);
  }
  sp->seq[sp->len] = 0;
  buf = strstr(sp->header, " LEN=");
  if (buf) {
    char *s;
    if (*(buf - 1) == ';') {
      buf -= 1;
      s = buf + 6;
      headerLen -= 6;
    } else {
      s = buf + 5;
      headerLen -= 5;
    }
    while (isdigit(*s)) {
      s += 1;
      headerLen -= 1;
    }
    while (*s)
      *buf++ = *s++;
  }
  buf = sp->header + headerLen - 1;
  while (iscntrl(*buf) || isspace(*buf))
    buf -= 1;
  res = snprintf(buf + 1, lenStr, "; LEN=%u\n", sp->len);
  if (res < 0 || res >= lenStr)
    fatal("Sequence too long: %u\n", sp->len);
  gc = ctr['G'] + ctr['C'];
  atgc = gc + ctr['A'] + ctr['T'];
  sp->GC_pct = (atgc == 0) ? 0.0 : 100.0 * (double) gc / (double) atgc;
  return 0;
}

static void
free_seq(seq_p_t sp)
{
  free(sp->seq);
  free(sp->header);
  free_buf(&sp->rb);
  if (sp->fName != NULL)
    close(sp->fd);
}

static void
seq_revcomp_inplace(seq_p_t seq)
{
  unsigned char *s = seq->seq;
  unsigned char *t = seq->seq + seq->len;
  unsigned char c;
  while (s < t) {
    c = dna_complement[*--t];
    *t = dna_complement[*s];
    *s++ = c;
  }
}

static void
init_col(col_p_t c, unsigned int size)
{
  c->size = size;
  c->nb = 0;
  if (size > 0)
    c->e.elt = (void **) xmalloc(size * sizeof(void *));
  else
    c->e.elt = NULL;
}

static void
add_col_elt(col_p_t c, void *elt, unsigned int grow)
{
  if (c->size <= c->nb) {
    c->size += grow;
    c->e.elt = (void **) xrealloc(c->e.elt, c->size * sizeof(void *));
  }
  c->e.elt[c->nb++] = elt;
}

#if 0 /* CI */
static void
add_unique_col_elt(col_p_t c, void *elt, unsigned int grow)
{
  unsigned int i;
  for (i = 0; i < c->nb; i++)
    if (c->e.elt[i] == elt)
      return;
  add_col_elt(c, elt, grow);
}

static void
merge_col(col_p_t c1, col_p_t c2)
{
  unsigned int i;
  for (i = 0; i < c2->nb; i++)
    add_col_elt(c1, c2->e.elt[i], COL_G_GROW);
}
#endif /* 0 CI */

static void
free_col(col_p_t c)
{
#ifndef NDEBUG
  memset(c->e.elt, 0, c->size * sizeof(void *));
#endif
  free(c->e.elt);
#ifndef NDEBUG
  memset(c, 0, sizeof(col_t));
#endif
}

static void
CreateMatrix(matrix_p_t m, signed char *data, unsigned int nElt)
{
  int i;
  unsigned int frame;
  unsigned int sSize = 1;
  unsigned int *step = (unsigned int *) xmalloc(sizeof(unsigned int)
						* m->order);
  unsigned int *sStep = (unsigned int *) xmalloc(sizeof(unsigned int)
						 * m->order);

  if (m->order < 1)
    fatal("CreateMatrix: order should be >=1 (%d)", m->order);
  m->m = (signed char **) xmalloc(sizeof(signed char *) * m->frames);
  /* Compute some stepping info.  */
  step[m->order - 1] = 4;
  sStep[m->order - 1] = 5;
  for (i = m->order - 2; i >= 0; i--) {
    step[i] = step[i + 1] * 4;
    sStep[i] = sStep[i + 1] * 5;
  }
  /* Check the size of the array.  */
  if (step[0] * m->frames != nElt)
    fatal("CreateMatrix: bad array size (%d, should be %d)",
	  nElt, step[0] * m->frames);
  /* Compute the score table size.  */
  for (frame = 0; frame < m->order; frame++)
    sSize *= 5;
  for (frame = 0; frame < m->frames; frame++) {
    signed char *ptr;
    /* Get space for the score tables.  */
    m->m[frame] = (signed char *) xmalloc(sizeof(signed char) * sSize);
    ptr = m->m[frame];
    /* Process the array.  */
    for (i = 0; i < (int) step[0]; i++) {
      int val = *data++;
      int j;
      /* Do not go below min.  */
      val = (val < options.min) ? options.min : val;
      *ptr++ = val;
      for (j = m->order - 1; j >= 0; j--) {
	if ((i + 1) % step[j] == 0) {
	  /* We have to fill in the next sStep[j]/5 score slots.  */
	  int k;
	  int stepping = sStep[j] / 5;
	  if (options.Nvalue == 0) {
	    /* Plain average thing.  */
	    for (k = 0; k < stepping - 1; k++) {
	      int avg = *(ptr - stepping) + *(ptr - stepping * 2)
		      + *(ptr - stepping * 3) + *(ptr - stepping * 4);
	      avg /= 4;
	      *ptr++ = avg;
	    }
	    *ptr++ = 0; /* Null expectation to accept an N.  */
	  } else {
	    /* Something a bit more funky...  */
	    for (k = 0; k < stepping - 1; k++) {
	      int avg;
	      int sorted[4];
	      sorted[0] = *(ptr - stepping);
	      sorted[1] = *(ptr - stepping * 2);
	      sorted[2] = *(ptr - stepping * 3);
	      sorted[3] = *(ptr - stepping * 4);
	      /* Need to sort the darn thing...  */
	      qsort(sorted, 4, sizeof(int), intCompare);
	      switch(options.Nvalue) {
	      case 1:
		avg = sorted[3];
		break;
	      case 2:
		avg = (sorted[3] + sorted[2]) / 2;
		break;
	      case 3:
		avg = (sorted[3] + sorted[2] + sorted[1]) / 3;
		break;
	      case -1:
		avg = sorted[0];
		break;
	      case -2:
		avg = (sorted[0] + sorted[1]) / 2;
		break;
	      case -3:
		avg = (sorted[0] + sorted[1] + sorted[2]) / 3;
		break;
	      default:
		fatal("Bad method (%d) to compute N score value.",
		      options.Nvalue);
	      }
	      *ptr++ = avg;
	    }
	    *ptr++ = 0; /* Null expectation to accept an N.  */
	  }
	}
      }
    }
  }
  free(step);
  free(sStep);
}

static inline unsigned int
GetCode(unsigned char c)
{
  switch (c) {
  case 'A':
    return 0; /* This is A.  */
  case 'C':
    return 1; /* This is C.  */
  case 'G':
    return 2; /* This is G.  */
  case 'T':
    return 3; /* This is T.  */
  }
  return 4; /* Everything else.  */
}

static inline void
findMax(int prev, int *prevV, int transit, int *bPrev, int *bScore)
{
  int score = prevV[prev] + transit;
  if (score > *bScore) {
    *bScore = score;
    *bPrev = prev;
  }
}

static inline void
findMax0(int prev, int *prevV, int *bPrev, int *bScore)
{
  int score = prevV[prev];
  if (score > *bScore) {
    *bScore = score;
    *bPrev = prev;
  }
}

static void
initIndices(int tsize, int startlen, int stoplen)
{
  unsigned int f;
  iBegin = -1;
  i5utr = 0;
  iStart = 1;
  iCds = iStart + startlen;
  iStop = iCds + 3;
  i3utr = iStop + stoplen;
  for (f = 0; f < 3; f++) {
    iInsAfter[f] = i3utr + f * tsize + 1;
    iDelAfter[f] = i3utr + (f + 3) * tsize - f + 1;
  }
  for (f = 0; f < 3; f++) {
    unsigned int f1;
    for (f1 = 0; f1 < 3; f1++) {
      if ((f1 + tsize)     % 3 == f)
	iInsNext[f] = iInsAfter[f1] + tsize - 1;
      if ((f1 + tsize + 1) % 3 == f)
	iDelNext[f] = iDelAfter[f1] + tsize - 2;
    }
  }
}

#ifdef DEBUG
static void
printIndex(unsigned int index, unsigned int len)
{
  char *s = (char *) xmalloc(sizeof(char) * (len + 1));
  int i;
  s[len] = 0;
  for (i = len - 1; i >= 0; i--) {
    int c = index % 5;
    index /= 5;
    switch (c) {
    case 0:
      s[i] = 'A';
      break;
    case 1:
      s[i] = 'C';
      break;
    case 2:
      s[i] = 'G';
      break;
    case 3:
      s[i] = 'T';
      break;
    default:
      s[i] = 'N';
    }
  }
  fputs(s, stderr);
  free(s);
}

static void
printInitStatus(unsigned int states, unsigned int seqLen,
		unsigned int tsize, unsigned int tableSize,
		unsigned int tindex,
		unsigned int *insTindex, unsigned int *delTindex)
{
  unsigned int f, i;
  fprintf(stderr, "Begin: %d\n", iBegin);
  fprintf(stderr, "5'UTR: %d\n", i5utr);
  fprintf(stderr, "Start: %d\n", iStart);
  fprintf(stderr, "Stop:  %d\n", iStop);
  fprintf(stderr, "3'UTR: %d\n", i3utr);
  for (f = 0; f < 3; f++)
    fprintf(stderr, "Frame %u: CDS %d, insert after/next %d/%d, delete after/next %d/%d\n",
	    f, iCds + f, iInsAfter[f], iInsNext[f], iDelAfter[f], iDelNext[f]);
  fprintf(stderr, "states %u, seq length %u, tsize %u, tableSize %u, tindex ",
	  states, seqLen, tsize, tableSize);
  printIndex(tindex, tsize);
  fprintf(stderr, "\ninsertion indices: ");
  for (i = 0; i < tsize; i++) {
    printIndex(insTindex[i], tsize);
    fprintf(stderr, " ");
  }
  fprintf(stderr, "\ndeletion indices: ");
  for (i = 0; i < tsize - 1; i++) {
    printIndex(delTindex[i], tsize);
    fprintf(stderr, " ");
  }
  fprintf(stderr, "\n");
}

static void
printCurrentStatus(unsigned int p, unsigned char c, unsigned int code,
		   unsigned int tindex, unsigned int tsize,
		   unsigned int *insTindex, unsigned int *delTindex,
		   unsigned int states,
		   int *currV, int *currTr)
{
  unsigned int i;
  fprintf(stderr, "%u:%c-%u: ", p, c, code);
  printIndex(tindex, tsize);
  fprintf (stderr, " /");
  for (i = 0; i < tsize; i++) {
    fprintf(stderr, " ");
    printIndex(insTindex[i], tsize);
  }
  fprintf(stderr, " /");
  for (i = 0; i < tsize - 1; i++) {
    fprintf(stderr, " ");
    printIndex(delTindex[i], tsize);
  }
  fprintf(stderr, "\n");
  for (i = 0; i < states; i++) {
    if (currV[i] < INT_MIN / 3)
      fprintf(stderr, " -inf/%2d", currTr[i]);
    else
      fprintf(stderr, "%5d/%2d", currV[i], currTr[i]);
    if ((i % 10) == 9)
      fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
}
#endif

/* returns frame if <state> is coding, relies on startoffset,
   startlength, stopoffset and stoplength to be multiples of
   3, returns -1 if not coding */
static inline int
getFrame(int state, int tsize, int startoffset, int stopoffset)
{
  int f;
  int d = state - iStart - startoffset + 1;
  if (d < 0)
    return -1;
  if (state < (iStop + stopoffset))
    return(d % 3);
  for (f = 1; f < tsize; f++) {
    if ((iInsAfter[(12 - f) % 3] + f) == state)
      return 0;
    if ((iInsAfter[(13 - f) % 3] + f) == state)
      return 1;
    if ((iInsAfter[(14 - f) % 3] + f) == state)
      return 2;
    if ((iDelAfter[(14 - f) % 3] + f - 1) == state)
      return 0;
    if ((iDelAfter[(15 - f) % 3] + f - 1) == state)
      return 1;
    if ((iDelAfter[(16 - f) % 3] + f - 1) == state)
      return 2;
  }
  return -1;
}

static int
Compute(seq_p_t seq, col_p_t mc, col_p_t rc, int reverse, int maxScore)
{
  matrix_p_t M[MT_COUNT];
  unsigned int i, code, f, s;
  unsigned char *p;
  int iCurr, bPrev, bScore, tmpPrev, tmpScore;
  unsigned int tableSize = 1;
  unsigned int tindex;
  unsigned int *insTindex, *delTindex;
  unsigned int states, mSize;
  int *currV, *prevV, *currTr;

  /* Find the right matrices.  */
  memset(M, 0, sizeof(M));
  for (i = 0; i < mc->nb; i ++) {
    if (seq->GC_pct >= mc->e.m[i]->CGmin
	&& seq->GC_pct <= mc->e.m[i]->CGmax
	&& M[mc->e.m[i]->matType] == NULL)
      M[mc->e.m[i]->matType] = mc->e.m[i];
  }
  for (i = 0; i < MT_COUNT; i++)
    if (M[i] == NULL)
      fatal("We have no %d matrix for %.2f GC in:\n %s",
	    i, seq->GC_pct, seq->header);
  /* initialize some more parameters */
  insTindex = (unsigned int *) xmalloc(sizeof(unsigned int)
				       * M[MT_CODING]->order);
  delTindex = (unsigned int *) xmalloc(sizeof(unsigned int)
				       * (M[MT_CODING]->order - 1));
  /* allocate tables and compute the state indices  */
  states = 2 + M[MT_START]->frames
	  + M[MT_STOP]->frames + 6 * M[MT_CODING]->order;
  mSize = sizeof(int) * seq->len * states;
  if (maxSize < mSize) {
    maxSize = mSize;
    V  = (int *) xrealloc(V,  maxSize);
    tr = (int *) xrealloc(tr, maxSize);
  }
  currV  = V; currTr = tr;
  /* size of score tables per frame */
  for (i = 0; i < M[MT_CODING]->order; i++)
    tableSize *= 5;
  /* tindex will point to the position representing the last
   * M[MT_CODING]->order chars on seq */
  /* tindex now represents all N's */
  tindex = tableSize - 1;
  for (i = 0; i < M[MT_CODING]->order - 1; i++)
    insTindex[i] = delTindex[i] = tindex;
  insTindex[i] = tindex;
  initIndices(M[MT_CODING]->order, M[MT_START]->frames, M[MT_STOP]->frames);
#ifdef DEBUG
  printInitStatus(states, seq->len, M[MT_CODING]->order, tableSize,
		  tindex, insTindex, delTindex);
#endif
  /* fill in the Viterbi and traceback tables, initialize for first char on seq */
  code = GetCode(seq->seq[0]);
  tindex = (5 * tindex + code) % tableSize;
  currV[i5utr] = options.ts5uPen + M[MT_UNTRANSLATED]->m[0][tindex];
  for (f = 0; f < M[MT_START]->frames; f++)
    currV[iStart+f] = options.min + M[MT_START]->m[f][code];
  for (f = 0; f <  3; f++)
    currV[iCds + f] = options.tscPen + M[MT_CODING]->m[f][tindex];
  for (f = 0; f < M[MT_STOP]->frames; f++)
    currV[iStop + f] = options.min + M[MT_STOP]->m[f][code];
  currV[i3utr] = options.ts3uPen + M[MT_UNTRANSLATED]->m[0][tindex];
  for (s = i3utr + 1; s < states; s++)
    currV[s] = INT_MIN / 2;
  for (s = 0; s < states; s++)
    currTr[s] = iBegin;
#ifdef DEBUG
  printCurrentStatus(0, seq->seq[0], code, tindex, M[MT_CODING]->order,
		     insTindex, delTindex, states, V, tr);
#endif
  /* fill in the Viterbi and traceback tables, main part */
  for (p = seq->seq + 1; *p; p++) {
    prevV = currV;
    currV += states;
    currTr += states;
    /* update index variables */
    code = GetCode(*p);
    for (i = M[MT_CODING]->order - 1; i > 0; i--)
      insTindex[i] = (5 * insTindex[i - 1] + code) % tableSize;
    if (M[MT_CODING]->order > 2)
      for (i = M[MT_CODING]->order - 2; i > 0; i--)
	delTindex[i] = (5 * delTindex[i - 1] + code) % tableSize;
    insTindex[0] = tindex;
    delTindex[0] = (25 * tindex + 20 + code) % tableSize;
    tindex = (5 * tindex + code) % tableSize;
    /* consider current nucleotide in 5'UTR */
    /* transitions UTR->UTR and CDS->CDS are presumed zero */
    currV[i5utr]  = prevV[i5utr] + M[MT_UNTRANSLATED]->m[0][tindex];
    currTr[i5utr] = i5utr;
    /* consider current nucleotide in start profile */
    currV[iStart]  = prevV[i5utr] + options.t5ucPen + M[MT_START]->m[0][code];
    currTr[iStart] = i5utr;
    iCurr = iStart;
    bPrev = iStart - 1;
    for (f = 1; f < M[MT_START]->frames; f++) {
      iCurr += 1;
      bPrev += 1;
      currV[iCurr] = prevV[bPrev] + M[MT_START]->m[f][code];
      currTr[iCurr] = bPrev;
    }
    /* consider current nucleotide in CDS */
    iCurr = iCds;
    bPrev = iCds - 1;
    bScore = prevV[bPrev];
    findMax(i5utr, prevV, options.min, &bPrev, &bScore);
    findMax0(iCurr + 2, prevV, &bPrev, &bScore);
    findMax0(iInsNext[0], prevV, &bPrev, &bScore);
    findMax0(iDelNext[0], prevV, &bPrev, &bScore);
    currV[iCurr] = bScore + M[MT_CODING]->m[0][tindex];
    currTr[iCurr] = bPrev;
    for (f = 1; f < 3; f++) {
      iCurr += 1;
      bPrev = iCurr - 1;
      bScore = prevV[bPrev];
      findMax0(iInsNext[f], prevV, &bPrev, &bScore);
      findMax0(iDelNext[f], prevV, &bPrev, &bScore);
      currV[iCurr] = bScore + M[MT_CODING]->m[f][tindex];
      currTr[iCurr] = bPrev;
    }
    /* consider current nucleotide in stop profile */
    bPrev = INT_MIN;
    bScore = INT_MIN;
    for (f = M[MT_START]->offset + 2; f < M[MT_START]->frames; f += 3)
      findMax(iStart + f, prevV, options.min, &bPrev, &bScore);
    for (f = 0; f < M[MT_CODING]->order; f++)
      findMax(iInsAfter[(14 - f) % 3] + f, prevV, options.min, &bPrev, &bScore);
    for (f = 0; f < M[MT_CODING]->order - 1; f++)
      findMax(iDelAfter[(15 - f) % 3] + f, prevV, options.min, &bPrev, &bScore);
    tmpPrev = bPrev;
    tmpScore = bScore;
    findMax(iCds + 2, prevV, options.tc3uPen, &bPrev, &bScore);
    currV[iStop] = bScore + M[MT_STOP]->m[0][code];
    currTr[iStop] = bPrev;
    iCurr = iStop;
    bPrev = iStop - 1;
    for (f = 1; f < M[MT_STOP]->frames; f++) {
      iCurr += 1;
      bPrev += 1;
      currV[iCurr] = prevV[bPrev] + M[MT_STOP]->m[f][code];
      currTr[iCurr] = bPrev;
    }
    /* consider current nucleotide in 3' UTR */
    bPrev = INT_MIN;
    bScore = INT_MIN;
    findMax0(i3utr - 1, prevV, &bPrev, &bScore);
    findMax0(i3utr, prevV, &bPrev, &bScore);
    findMax(iCds+2, prevV, options.min, &bPrev, &bScore);
    currV[i3utr]  = bScore + M[MT_UNTRANSLATED]->m[0][tindex];
    currTr[i3utr] = bPrev;
    /* consider current nucleotide in CDS after insertion */
    for (f = 0; f < 3; f++) {
      iCurr = iInsAfter[f];
      bPrev = iCds + f;
      currV[iCurr] = prevV[bPrev] + options.iPen;
      currTr[iCurr] = bPrev;
      iCurr = iInsAfter[f];
      bPrev = iCurr - 1;
      for (i = 1; i < M[MT_CODING]->order; i++) {
	iCurr += 1;
	bPrev += 1;
	currV[iCurr] = prevV[bPrev] + M[MT_CODING]->m[(i + f) % 3][insTindex[i]];
	currTr[iCurr] = bPrev;
      }
    }
    /* consider current nucleotide in CDS after deletion */
    for (f = 0; f < 3; f++) {
      iCurr = iDelAfter[f];
      bPrev = iCds + f;
      currV[iCurr] = prevV[bPrev] + options.dPen
		    + M[MT_CODING]->m[(f + 2) % 3][delTindex[0]];
      currTr[iCurr] = bPrev;
      iCurr = iDelAfter[f];
      bPrev = iCurr - 1;
      for (i = 1; i < M[MT_CODING]->order - 1; i++) {
	iCurr += 1;
	bPrev += 1;
	currV[iCurr] = prevV[bPrev] + M[MT_CODING]->m[(i + f + 2) % 3][delTindex[i]];
	currTr[iCurr] = bPrev;
      }
    }
#ifdef DEBUG
    printCurrentStatus(p - seq->seq, *p, code, tindex, M[MT_CODING]->order,
		       insTindex, delTindex, states, currV, currTr);
#endif
  }
  /* fill in the Viterbi and traceback tables, terminate and find best */
  bPrev = i5utr;
  bScore = currV[bPrev] + options.t5uePen;
  for (f = 0; f < M[MT_START]->frames; f++)
    findMax(iStart+f, currV, options.min, &bPrev, &bScore);
  for (f = 0; f < 3; f++) {
    findMax(iCds+f, currV, options.tcePen, &bPrev, &bScore);
    for (i = 0; i < M[MT_CODING]->order; i++)
      findMax(iInsAfter[f] + i, currV, options.tcePen, &bPrev, &bScore);
    for (i = 0; i < M[MT_CODING]->order - 1; i++)
      findMax(iDelAfter[f] + i, currV, options.tcePen, &bPrev, &bScore);
  }
  for (f = 0; f < M[MT_STOP]->frames; f++)
    findMax(iStop+f, currV, options.min, &bPrev, &bScore);
  findMax(i3utr, currV, options.t3uePen, &bPrev, &bScore);
#ifdef DEBUG
  fprintf(stderr, "finished to fill Viterbi matrix, best score %d in state %d\n",
	  bScore, bPrev);
#endif
  /* traceback and generate coding sequences starting from bPrev (confidence bScore) */
  iCurr = bPrev;
  p -= 1;
  while(iCurr != iBegin) {
    int iOld = -1, rStart, rStop;
    unsigned char *r, *q;
    /* skip non coding */
    while((iBegin < iCurr
	   && iCurr < iStart + M[MT_START]->offset - 1)
	  || (iStop + M[MT_STOP]->offset - 1 < iCurr
	      && iCurr < iInsAfter[0])) {
#ifdef DEBUG
      fprintf(stderr, "trace back non-coding: state %d position %4d(%c)\n",
	      iCurr, p - seq->seq, *p);
#endif
      iCurr = currTr[iCurr];
      currTr -= states;
      p -= 1;
    }
    /* handle coding */
    if (iCurr != iBegin) {
      unsigned char *res = (unsigned char *) xmalloc(sizeof(unsigned char)
						     * 2 * seq->len);
      int rScore = V[(p - seq->seq) * states + iCurr];
      r = res;
      rStop = (p - seq->seq);
      if (getFrame(iCurr, M[MT_CODING]->order,
		   M[MT_START]->offset, M[MT_STOP]->offset) == 0) {
	*r++ = 'X';
	*r++ = 'X';
      }
      if (getFrame(iCurr, M[MT_CODING]->order,
		   M[MT_START]->offset, M[MT_STOP]->offset) == 1)
	*r++ = 'X';
      while((iStart + M[MT_START]->offset - 1 <= iCurr
	     && iCurr <= iStop + M[MT_STOP]->offset - 1)
	    || iInsAfter[0] <= iCurr) {
	int done = 0;
	for (f = 0; f < 3; f++) {
	  if (iCurr == iInsAfter[f]) {
	    *r++ = tolower(*p);
	    done = 1;
	  }
	  if (iCurr == iDelAfter[f]) {
	    *r++ = toupper(*p);
	    *r++ = 'X';
	    done = 1;
	  }
	}
	if (done == 0)
	  *r++ = toupper(*p);
	/* remove stop-profile penalty from coding score */
	if (iCurr == iCds + 2 && iOld == iStop)
	  rScore -= options.tc3uPen;
#ifdef DEBUG
      fprintf(stderr, "trace back     coding: state %2d(%2d) position %4d(%c)\n",
	      iCurr, getFrame(iCurr, M[MT_CODING]->order,
	      M[MT_START]->offset, M[MT_STOP]->offset), p-seq->seq, *(r-1));
#endif
	iOld = iCurr;
	iCurr = currTr[iCurr];
	currTr -= states;
	p -= 1;
      }
      rStart = p - seq->seq + 1;
      if (p >= seq->seq)
	rScore -= V[(p - seq->seq) * states + iCurr];
      if (rScore > maxScore)
	maxScore = rScore;
      if (getFrame(iOld, M[MT_CODING]->order,
		   M[MT_START]->offset, M[MT_STOP]->offset) == 1)
	*r++ = 'X';
      if (getFrame(iOld, M[MT_CODING]->order,
		   M[MT_START]->offset, M[MT_STOP]->offset) == 2) {
	*r++='X';
	*r++='X';
      }
      *r-- = 0;
      /* reverse the array and add to the result-array */
      q = res;
      while (q < r) {
	unsigned char c = *r;
	*r-- = *q;
	*q++ = c;
      }
#ifdef DEBUG
      fprintf(stderr, "found coding %s, add to results, state %d \n",
	      res, iCurr);
#endif
      if (rStop - rStart >= options.minLen) {
	result_p_t r = (result_p_t) xmalloc(sizeof(result_t));
	r->score = rScore;
	r->start = rStart;
	r->stop = rStop;
	r->reverse = reverse;
	r->s = res;
	add_col_elt(rc, r, 8);
      } else
	free(res);
    }
  }
  /* clean up */
  free(insTindex);
  free(delTindex);
  return maxScore;
}

static void
LoadMatrix(const char *fName, col_p_t mc)
{
  read_buf_t rb;
  char *buf;
  int fd = open(fName, O_RDONLY);
  if (fd == -1)
    fatal("Could not open file %s: %s(%d)\n", fName, strerror(errno), errno);
  init_col(mc, 16);
  init_buf(&rb);
  buf = read_line_buf(&rb, fd);
  while (rb.lc > 0) {
    if (strncmp(buf, "FORMAT: ", 8) == 0) {
      matrix_p_t m = (matrix_p_t) xmalloc(sizeof(matrix_t));
      char name[256], fType[256], mType[256];
      int res;
      unsigned int size = 256;
      signed char *data = (signed char *) xmalloc(sizeof(signed char) * size);
      unsigned int nElt = 0;
      res = sscanf(buf,
		   "FORMAT: %255s %255s %255s %u %u %d s C+G: %lf %lf",
		   name, fType, mType, &m->order, &m->frames, &m->offset,
		   &m->CGmin, &m->CGmax);
      if (m->CGmin < 0.0)
	m->CGmin = 0.0;
      if (m->CGmin > 0.0)
	m->CGmin += options.percent;
      m->CGmax += options.percent;
      if (m->CGmax > 100.0)
	m->CGmax = 100.0;
      if (res != 8 || (buf = read_line_buf(&rb, fd)) == NULL || rb.lc == 0)
	fatal("Bad data header format in file %s, near %s (%d)\n",
	      fName, name, res);
      m->name = strdup(name);
      m->kind = strdup(mType);
      m->matType = MT_UNKNOWN;
      if (strncmp(fType, "CODING", 6) == 0)
	m->matType = MT_CODING;
      if (strncmp(fType, "UNTRANSLATED", 12) == 0)
	m->matType = MT_UNTRANSLATED;
      if (strncmp(fType, "START", 5) == 0)
	m->matType = MT_START;
      if (strncmp(fType, "STOP", 4) == 0)
	m->matType = MT_STOP;
      while (buf[0] == '-' || isdigit(buf[0])) {
	int a, c, g, t;
	res = sscanf(buf, "%d %d %d %d", &a, &c, &g, &t);
	if (res != 4 || (buf = read_line_buf(&rb, fd)) == NULL)
	  fatal("Bad data format in file %s, near %s (%d)\n",
		fName, name, res);
	if (nElt + 4 > size) {
	  size += 256;
	  data = (signed char *) xrealloc(data, sizeof(signed char) * size);
	}
	data[nElt++] = a;
	data[nElt++] = c;
	data[nElt++] = g;
	data[nElt++] = t;
      }
      CreateMatrix(m, data, nElt);
      add_col_elt(mc, m, 16);
      free(data);
    } else if ((buf = read_line_buf(&rb, fd)) == NULL)
      fatal("Probable bug in read_line_buf while reading %s\n", fName);
  }
  close(fd);
  free_buf(&rb);
}

static void
remove_lc(unsigned char *s)
{
  unsigned char *t = s;
  while (*s) {
    if (isupper(*s))
      *t++ = *s;
    s += 1;
  }
  *t = 0;
}

static char *
na2aa(unsigned char *s)
{
  static char *CABC = "KNKNTTTTRSRSIIMI"  /* AAA AAC ... ATT */
		      "QHQHPPPPRRRRLLLL"  /* CAA CAC ... CTT */
		      "EDEDAAAAGGGGVVVV"  /* GAA GAC ... GTT */
		      "OYOYSSSSOCWCLFLF"; /* TAA TAC ... TTT */
  static char *CNBC = "XTXXXPRLXAGVXSXX"; /* AAN ACN ... TTN */
  char *res = xmalloc(sizeof(char) * (strlen((char *) s) / 3 + 2));
  char *cur = res;

  while (*s) {
    int idx = 0;
    /* Check first nt.  */
    switch (*s) {
    case 'A':
      break;
    case 'C':
      idx = 1;
      break;
    case 'G':
      idx = 2;
      break;
    case 'T':
      idx = 3;
      break;
    default:
      idx = -1;
    }
    s += 1;
    if (*s == 0) {
      *cur++ = 'X';
      break;
    }
    if (idx == -1) {
      *cur++ = 'X';
      s += 1;
      if (*s == 0)
	break;
      s += 1;
      continue;
    }
    idx <<= 2;
    /* Check second nt.  */
    switch (*s) {
    case 'A':
      break;
    case 'C':
      idx += 1;
      break;
    case 'G':
      idx += 2;
      break;
    case 'T':
      idx += 3;
      break;
    default:
      idx = -1;
    }
    s += 1;
    if (*s == 0) {
      if (idx == -1)
	*cur++ = 'X';
      else
	*cur++ = CNBC[idx];
      break;
    }
    if (idx == -1) {
      *cur++ = 'X';
      s += 1;
      continue;
    }
    idx <<= 2;
    /* Check third nt.  */
    switch (*s) {
    case 'A':
      break;
    case 'C':
      idx += 1;
      break;
    case 'G':
      idx += 2;
      break;
    case 'T':
      idx += 3;
      break;
    default:
      *cur++ = CNBC[idx >> 2];
      idx = -1;
    }
    if (idx != -1)
      *cur++ = CABC[idx];
    s += 1;
  }
  *cur = 0;
  return res;
}

static void
showResults(col_p_t rc, seq_p_t seq, unsigned char *rSeq, int maxScore)
{
  unsigned int i;
  unsigned int cnt = 0;
  if (options.maxOnly != 0) {
    result_p_t r = NULL;
    if (options.out == NULL)
      return;
    for (i = 0; i < rc->nb; i++)
      if (rc->e.r[i]->score == maxScore) {
	r = rc->e.r[i];
	break;
      }
    i = 0;
    while (!isspace(seq->header[i]))
      i += 1;
    if (r != NULL)
      fprintf(options.out, "%.*s %d %u %u %u %c\n", i, seq->header,
	      r->score, r->start + 1, r->stop + 1, seq->len,
	      r->reverse ? '-' : '+');
    else
      fprintf(options.out, "%.*s %d\n", i, seq->header, maxScore);
    return;
  }
  if (options.all != 0) {
    fprintf(stderr, "Sorry, options.all is unimplemented yet...\n");
/*	my $scores = "";
	my $outSeq = "";
	my $lastPos = 0;
	foreach my $r (@$rres) {
	    my $curScore = $$r[0];
	    $scores .= $curScore . " ";
	    if ($theRealMax * $main::both > $curScore) { next; }
	    if ($$r[1] > $lastPos) {
		my $s = substr($seq->{_seq}, $lastPos, $$r[1] - $lastPos);
		$s =~ tr/A-Z/a-z/;
		$outSeq .= $s;
	    }
	    $outSeq .= $$r[3];
	    $lastPos = $$r[2] + 1;
	}
	if ($lastPos < $seq->seqLength) {
	    my $s = substr($seq->{_seq}, $lastPos);
	    $s =~ tr/A-Z/a-z/;
	    $outSeq .= $s;
	}
	my $head = $seq->seqHead;
	$head =~ s/^(\S+)/$1 $scores/;
	print $main::out "$head\n";
	$outSeq =~ s/(.{$main::sWidth})/$1\n/g;
	$outSeq =~ s/\s+$//; # remove a trailing newline, since we add one below.
	print $main::out "$outSeq\n"; */
    return;
  }
  for (i = 0; i < rc->nb; i++) {
    result_p_t r = rc->e.r[i];
    char *h = seq->header;
    size_t len = strlen(h);
    char *buf, *ptr;
    if ((double) maxScore * options.both > (double) r->score)
      continue;
    buf = (char *) xmalloc((len + 256) * sizeof(char));
    ptr = buf;
    while (*h && *h != '|' && !isspace(*h))
      *ptr++ = *h++;
    if (*h == '|') {
      while (*h && *h != '|' && !isspace(*h))
	*ptr++ = *h++;
    }
    if (cnt > 0)
      *ptr++ = 'a' + cnt - 1;
    cnt += 1;
    while (*h && !isspace(*h))
      *ptr++ = *h++;
    sprintf(ptr, " %d %d %d %s", r->score, r->start + 1, r->stop + 1, h);
    if (r->reverse != 0) {
      char *t = strstr(buf, "; minus strand");
      if (t != NULL) {
	while (*(t + 14) != 0) {
	  *t = *(t + 14);
	  t += 1;
	}
	*t = 0;
      } else {
	len = strlen(buf);
	while (isspace(buf[len - 1]))
	  len -= 1;
	strcpy(buf + len, "; minus strand\n");
      }
    }
    if (options.transl != NULL) {
      char *ps;
      remove_lc(r->s);
      ps = na2aa(r->s);
      len = strlen(buf);
      while (isspace(buf[len - 1]))
	len -= 1;
      fprintf(options.transl, "%.*s; translated\n", (unsigned int) len, buf);
      len = strlen(ps);
      ptr = ps;
      /* remove trailing stop codon(s).  */
      while (len > 0 && ps[len - 1] == 'O') {
	len -= 1;
	ps[len] = 0;
      }
      /* translate intermediate stop codons by X.  */
      while (*ptr) {
	if (*ptr == 'O')
	  *ptr = 'X';
	ptr += 1;
      }
      ptr = ps;
      while (len > options.sWidth) {
	fprintf(options.transl, "%.*s\n", options.sWidth, ptr);
	len -= options.sWidth;
	ptr += options.sWidth;
      }
      fprintf(options.transl, "%s\n", ptr);
      free(ps);
    }
    if (options.out != NULL) {
      fputs(buf, options.out);
      if (options.no_del != 0)
	remove_lc(r->s);
      len = strlen((char *) r->s);
      ptr = (char *) r->s;
      while (len > options.sWidth) {
	fprintf(options.out, "%.*s\n", options.sWidth, ptr);
	len -= options.sWidth;
	ptr += options.sWidth;
      }
      fprintf(options.out, "%s\n", ptr);
    }
    free(buf);
  }
}

static void
process_file(const char *fName, col_p_t mc)
{
  seq_t seq;
  col_t rc;
  init_col(&rc, 8);
  init_seq(fName, &seq);
  while (get_next_seq(&seq) == 0) {
    unsigned int i;
    int maxScore = INT_MIN;
    unsigned char *rSeq = NULL;
    if (seq.len == 0)
      continue;
    maxScore = Compute(&seq, mc, &rc, 0, maxScore);
    if (options.single == 0) {
      if (options.all)
	rSeq = (unsigned char *) strdup((char *) seq.seq);
      seq_revcomp_inplace(&seq);
      maxScore = Compute(&seq, mc, &rc, 1, maxScore);
      if (options.all) {
	unsigned char *tem = rSeq;
	rSeq = seq.seq;
	seq.seq = tem;
      }
    }
    showResults(&rc, &seq, rSeq, maxScore);
    for (i = 0; i < rc.nb; i++) {
      result_p_t r = rc.e.r[i];
      free(r->s);
      free(r);
    }
    rc.nb = 0;
    free(rSeq);
  }
  free_seq(&seq);
  free_col(&rc);
/*
    my $bigMax = ESTScan::Compute($seq->{_seq}, $main::iPen, $main::dPen, $main::min,
				  $main::maxOnly == 0 ? \@res : undef, $matIndex,
				  $utrMatIndex, $startMatIndex, $stopMatIndex,
				  $main::minLen);
    my $result_nr = -1;
    if ($main::single != 0) {
	showResults($seq, \@res, $bigMax, \$result_nr, $bigMax);
	next;
    }
    my @resRev;
    my $seqRev = $seq->revComp;
    my $bigMaxRev = ESTScan::Compute($seqRev->{_seq},$main::iPen,$main::dPen,$main::min,
				     $main::maxOnly == 0 ? \@resRev : undef, $matIndex,
				     $utrMatIndex, $startMatIndex, $stopMatIndex,
				     $main::minLen);

    my $theRealMax = $bigMax >= $bigMaxRev ? $bigMax : $bigMaxRev;
    showResults($seq, \@res, $bigMax, \$result_nr, $theRealMax);
    showResults($seqRev, \@resRev, $bigMaxRev, \$result_nr, $theRealMax);
}
*/
}

int
main(int argc, char *argv[])
{
  const char *ESTScanDir;
  int getHelp = 0;
  col_t mc;
#ifdef DEBUG
  unsigned int i;
  mcheck(NULL);
  mtrace();
#endif
  argv0 = argv[0];
  if (setlocale(LC_ALL, "POSIX") == NULL)
    fprintf(stderr, "%s: Warning: could not set locale to POSIX\n", argv[0]);
  /* Default options.  */
  ESTScanDir = getenv("ESTSCANDIR");
  if (ESTScanDir == NULL)
    ESTScanDir = "/usr/molbio/share/ESTScan";
  options.matrix = xmalloc((strlen(ESTScanDir) + 9) * sizeof(char));
  strcat(strcpy(options.matrix, ESTScanDir), "/Hs.smat");
  options.min = -100;
  options.dPen = -50;
  options.iPen = -50;
  options.ts5uPen = -10;
  options.tscPen = -10;
  options.ts3uPen = -5;
  options.t5ucPen = -80;
  options.t5uePen = -40;
  options.tc3uPen = -80;
  options.tcePen = -40;
  options.t3uePen = -20;
  options.percent = 4.0;
  options.Nvalue = 0;
  options.sWidth = 60;
  options.all = 0;
  options.maxOnly = 0;
  options.skipLen = 1;
  options.minLen = 50;
  options.transl = NULL;
  options.out = stdout;
  options.both = 1.0;
  options.no_del = 0;
  options.single = 0;
  while (1) {
    int c = getopt(argc, argv, "ab:d:hi:l:M:m:N:nOo:p:Ss:T:t:vw:");
    if (c == -1)
      break;
    switch (c) {
    case 'a':
      options.all = 1;
      break;
    case 'b':
      options.both = atof(optarg);
      break;
    case 'd':
      options.dPen = atoi(optarg);
      break;
    case 'h':
      getHelp = 1;
      break;
    case 'i':
      options.iPen = atoi(optarg);
      break;
    case 'l':
      options.minLen = atoi(optarg);
      break;
    case 'M':
      options.matrix = optarg;
      break;
    case 'm':
      options.min = atoi(optarg);
      break;
    case 'N':
      options.Nvalue = atoi(optarg);
      break;
    case 'n':
      options.no_del = 1;
      break;
    case 'O':
      options.maxOnly = 1;
      break;
    case 'o':
      if (strcmp(optarg, "-") == 0) {
	options.out = stdout;
      } else {
	options.out = fopen(optarg, "w");
	if (options.out == NULL)
	  fatal("Couldn't create file %s: %s (%d)\n", optarg,
		strerror(errno), errno);
      }
      break;
    case 'p':
      options.percent = atof(optarg);
      break;
    case 'S':
      options.single = 1;
      break;
    case 's':
      options.skipLen = atoi(optarg);
      break;
    case 'T':
/* if (defined $opt{'T'}) {
    my $nbProbs = ($opt{'T'} =~ s/,/,/g) + 1;
    if (($nbProbs) != 8) {
	usage("Wrong number of transition probabilities("  . $nbProbs . ")");
    }
    ($main::ts5uPen,$main::tscPen,$main::ts3uPen,$main::t5ucPen,
     $main::t5uePen,$main::tc3uPen,$main::tcePen,$main::t3uePen) = split(/,/,$opt{'T'});
} */
      fputs("Option T is not yet implemented...\n", stderr);
      return 1;
    case 't':
      if (strcmp(optarg, "-") == 0) {
	options.transl = stdout;
	if (options.out == options.transl)
	  options.out = NULL;
      } else {
	options.transl = fopen(optarg, "w");
	if (options.transl == NULL)
	  fatal("Couldn't create file %s: %s (%d)\n", optarg,
		strerror(errno), errno);
      }
      break;
    case 'v':
      fputs(Version, stderr);
      return 0;
    case 'w':
      options.sWidth = atoi(optarg);
      break;
    case '?':
      break;
    default:
      fprintf(stderr, "?? getopt returned character code 0%o ??\n", c);
    }
  }
  if (getHelp) {
    fprintf(stderr, Usage, argv[0], options.both, options.dPen,
	    options.iPen, options.minLen, options.matrix, options.min,
	    options.Nvalue, options.percent, options.skipLen,
	    options.ts5uPen, options.tscPen, options.ts3uPen,
	    options.t5ucPen, options.t5uePen, options.tc3uPen,
	    options.tcePen, options.t3uePen, options.sWidth);
    return 1;
  }
  LoadMatrix(options.matrix, &mc);
#ifdef DEBUG
  fprintf(stderr, "We have loaded %u matrices:\n", mc.nb);
  for (i = 0; i < mc.nb; i++) {
    matrix_p_t m = mc.e.m[i];
    fprintf(stderr, "%u: %s %s %d %.2f %.2f %u %u %d\n", i,
	    m->name, m->kind, m->matType, m->CGmin, m->CGmax,
	    m->order, m->frames, m->offset);
  }
#endif
  if (optind >= argc)
    process_file(NULL, &mc);
  else
    while (optind < argc)
      process_file(argv[optind++], &mc);
#ifdef DEBUG
  for (i = 0; i < mc.nb; i++) {
    matrix_p_t m = mc.e.m[i];
    unsigned int j;
    for (j = 0; j < m->frames; j++)
      free(m->m[j]);
    free(m->m);
    free(m->name);
    free(m->kind);
    free(m);
  }
  free_col(&mc);
  free(V);
  free(tr);
  free(options.matrix);
#endif
  return 0;
}
