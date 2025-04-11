/* $Id: estscan_.c,v 1.2 2006/12/14 23:59:58 c4chris Exp $
 *
 * Christian Iseli, LICR ITO, Christian.Iseli@licr.org
 *
 * Copyright (c) 1999 Swiss Institute of Bioinformatics. All rights reserved.
 */
#include "estscan_.h"
#include <stdio.h>
#include <limits.h>
#ifndef __GNUC__
#include <alloca.h>
#endif
#if !defined(__GNUC__) && defined(sun)
#define inline
#endif

#define __noDEBUG

static int
intCompare(const void *a, const void *b)
{
  return (*(int *)a < *(int *)b) ? -1 : ((*(int *)a > *(int *)b) ? 1 : 0);
}

/* At this point, we do not allow the deletion of unused matrices...  */

#define CODING 0
#define UNTRANSLATED 1
#define START 2
#define STOP 3

typedef struct _matrixGroup {
  signed char **m;
  int matType;
  int order;
  int frames;
  int offset;
  int CGmin;
  int CGmax;
} matrixGroup_t;

#define ALLOC_CHUNK 10

static matrixGroup_t *mg = NULL;
static int groupMax = 0;
static int nextGroup = 0;

static void
checkMGSpace(void)
{
  if (mg == 0) {
    /* First allocation.  */
    mg = (matrixGroup_t *) malloc(sizeof(matrixGroup_t) * ALLOC_CHUNK);
    groupMax += ALLOC_CHUNK;
  }
  if (nextGroup >= groupMax) {
    /* Grow the thing.  */
    mg = (matrixGroup_t *) realloc(mg, sizeof(matrixGroup_t) * (groupMax + ALLOC_CHUNK));
    groupMax += ALLOC_CHUNK;
  }
}

#ifdef __DEBUG
static void
printIndex(int index, int len)
{
  char *s = (char *)malloc(sizeof(char) * (len + 1));
  int i;
  s[len] = 0;
  for (i = len - 1; i >= 0; i--) {
    int c = index % 5;
    index /= 5;
    switch (c) {
    case 0:  s[i] = 'A'; break;
    case 1:  s[i] = 'C'; break;
    case 2:  s[i] = 'G'; break;
    case 3:  s[i] = 'T'; break;
    default: s[i] = 'N';
    }
  }
  printf("%s", s);
  free(s);
}
#endif /* __DEBUG */

int
CreateMatrix(int matType, int order, int frames, int offset, int CGmin, int CGmax,
	     SV *rv, int min, int method)
{
  AV *av;
  int avlen, i, frame;
  int sSize = 1;
  int *step = (int *)alloca(sizeof(int) * order);
  int *sStep = (int *)alloca(sizeof(int) * order);
  signed char **fScore = (signed char **)alloca(sizeof(signed char *) * frames);

  if ( SvROK( rv ) && (SvTYPE(SvRV(rv)) == SVt_PVAV) ) av = (AV*)SvRV(rv);
  else { warn("CreateMatrix: 2nd parameter was not an array ref"); return -1; }
  if (order < 1) { warn("CreateMatrix: order should be >=1 (%d)", order); return -1; }

  /* Compute some stepping info.  */
  step[order - 1] = 4;
  sStep[order - 1] = 5;
  for (i = order - 2; i >= 0; i--) {
    step[i] = step[i + 1] * 4;
    sStep[i] = sStep[i + 1] * 5;
  }
  /* The length of the array passed in parameter.  */
  avlen = av_len( av ) + 1;
  /* Check the size of the array.  */
  if (step[0] * frames != avlen) {
    warn("CreateMatrix: bad array size (%d, should be %d)", avlen, step[0] * frames);
    return -1;
  }
  /* Compute the score table size.  */
  for (i = 0; i < order; i++) sSize *= 5;
  for (frame = 0; frame < frames; frame++) {
    signed char *ptr;
    /* Get space for the score tables.  */
    fScore[frame] = (signed char *)malloc(sSize);
    ptr = fScore[frame];
    /* Process the array.  */
    for (i = 0; i < step[0]; i++) {
      SV *sv = av_shift(av);
      int val = SvIV(sv);
      int j;
      /* Do not go below min.  */
      val = (val < min) ? min : val;
      *ptr++ = val;
      for (j = order - 1; j >= 0; j--) {
	if ((i + 1) % step[j] == 0) {
	  /* We have to fill in the next sStep[j]/5 score slots.  */
	  int k;
	  int stepping = sStep[j] / 5;
	  if (method == 0) {
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
	      int avg = 0;
	      int sorted[4];
	      sorted[0] = *(ptr - stepping);
	      sorted[1] = *(ptr - stepping * 2);
	      sorted[2] = *(ptr - stepping * 3);
	      sorted[3] = *(ptr - stepping * 4);
	      /* Need to sort the darn thing...  */
	      qsort(sorted, 4, sizeof(int), intCompare);
	      switch(method) {
	      case 1: avg = sorted[3]; break;
	      case 2: avg = (sorted[3] + sorted[2]) / 2; break;
	      case 3: avg = (sorted[3] + sorted[2] + sorted[1]) / 3; break;
	      case -1: avg = sorted[0]; break;
	      case -2: avg = (sorted[0] + sorted[1]) / 2; break;
	      case -3: avg = (sorted[0] + sorted[1] + sorted[2]) / 3; break;
	      default: die("Bad method (%d) to compute N score value.", method);
	      }
	      *ptr++ = avg;
	    }
	    *ptr++ = 0; /* Null expectation to accept an N.  */
	  }
	}
      }
    }
  }
  /* Prepare room to store the matrix pointers.  */
  checkMGSpace();
  mg[nextGroup].m = malloc(sizeof(signed char *) * frames);
  for (frame = 0; frame < frames; frame++) {
    mg[nextGroup].m[frame] = fScore[frame];
  }
  mg[nextGroup].matType = matType;
  mg[nextGroup].order = order;
  mg[nextGroup].frames = frames;
  mg[nextGroup].offset = offset;
  mg[nextGroup].CGmin = CGmin;
  mg[nextGroup].CGmax = CGmax;  
  return nextGroup++;
}


/* transition log-probabilities */
static int ts5uPen;
static int tscPen;
static int ts3uPen;
static int t5ucPen;
static int t5uePen;
static int tc3uPen;
static int tcePen;
static int t3uePen;

int
StoreTransits(int ts5u, int tsc, int ts3u, int t5uc, 
	      int t5ue, int tc3u, int tce, int t3ue)
{
  ts5uPen = ts5u;
  tscPen = tsc;
  ts3uPen = ts3u;
  t5ucPen = t5uc;
  t5uePen = t5ue;
  tc3uPen = tc3u;
  tcePen = tce;
  t3uePen = t3ue;  
  return 0;
}

inline int
GetCode(unsigned char c)
{
  c &= 0x1f; /* Get lower bits, get rid of upper/lower info.  */
  switch (c) {
  case 1:  return 0; /* This is A.  */
  case 3:  return 1; /* This is C.  */
  case 7:  return 2; /* This is G.  */
  case 20: return 3; /* This is T.  */
  default: return 4; /* Everything else.  */
  }
}

double
ComputeGC(const char *s)
     /* determine matrix indices according to GC-content */
{
  int ctr[256], gc, atgc;
  ctr['a'] = 0; ctr['A'] = 0; ctr['c'] = 0; ctr['C'] = 0; 
  ctr['g'] = 0; ctr['G'] = 0; ctr['t'] = 0; ctr['T'] = 0; 
  while(*s) ctr[(unsigned char) *s++]++;
  gc = ctr['g'] + ctr['G'] + ctr['c'] + ctr['C']; 
  atgc = gc + ctr['a'] + ctr['A'] + ctr['t'] + ctr['T']; 
  return (atgc == 0) ? 0.0 : 100.0*(float)gc/(float)atgc;
}

inline void 
findMax(int prev, int *prevV, int transit, int *bPrev, int *bScore)
{ 
  int score = prevV[prev] + transit;
  if (score > *bScore) { *bScore = score; *bPrev = prev; } 
}

inline void 
findMax0(int prev, int *prevV, int *bPrev, int *bScore)
{ 
  int score = prevV[prev];
  if (score > *bScore) { *bScore = score; *bPrev = prev; } 
}

/* declaration of indexes also used in getFrame */
static int iBegin, i5utr, iStart, iCds, iStop, i3utr;
/* next tsize states implement insertion/deletion after nucleotide in frame index */
static int iInsAfter[3], iDelAfter[3];
/* last tsize states implemented insertion/deletion before nucleotide in frame index */
static int iInsNext[3], iDelNext[3];

static void 
initIndices(int tsize, int startlen, int stoplen)
{
  int f, f1;
  iBegin = -1; i5utr = 0; iStart = 1; iCds = iStart + startlen; 
  iStop = iCds + 3; i3utr = iStop + stoplen;
  for (f=0; f<3; f++) { 
    iInsAfter[f] = i3utr + f*tsize + 1; 
    iDelAfter[f] = i3utr + (f + 3)*tsize - f + 1; 
  }
  for (f = 0; f < 3; f++) { 
    for (f1 = 0; f1 < 3; f1++) { 
      if ((f1 + tsize)     % 3 == f) iInsNext[f] = iInsAfter[f1] + tsize - 1; 
      if ((f1 + tsize + 1) % 3 == f) iDelNext[f] = iDelAfter[f1] + tsize - 2; 
    }
  }
}

#ifdef __DEBUG
static void 
printInitStatus(int cdsIndex, int utrIndex, int startIndex, int stopIndex, 
		int states, int seqLen, int tsize, int tableSize, int tindex, 
		int *insTindex, int *delTindex) 
{
  int f, i;
  printf("Begin: %d\n", iBegin);
  printf("5'UTR: %d\n", i5utr);
  printf("Start: %d\n", iStart);
  printf("Stop:  %d\n", iStop);
  printf("3'UTR: %d\n", i3utr);
  for (f = 0; f < 3; f++) 
    printf("Frame %d: CDS %d, insert after/next %d/%d, delete after/next %d/%d\n",
	   f, iCds + f, iInsAfter[f], iInsNext[f], iDelAfter[f], iDelNext[f]);
  printf("matrix indices: %d, %d, %d, %d\n", cdsIndex, utrIndex, startIndex, stopIndex);
  printf("states %d, seq length %d, tsize %d, tableSize %d, tindex ", 
	 states, seqLen, tsize, tableSize);
  printIndex(tindex, tsize);
  printf("\ninsertion indices: ");
  for (i = 0; i < tsize; i++) { printIndex(insTindex[i], tsize); printf(" "); }
  printf("\ndeletion indices: ");
  for (i = 0; i < tsize - 1; i++) { printIndex(delTindex[i], tsize); printf(" "); }
  printf("\n");
}

static void
printCurrentStatus(int p, char c, int code, int tindex, int tsize, 
		   int *insTindex, int *delTindex, int states, int *currV, int *currTr)
{
  int i;
  printf("%d:%c-%d: ", p, c, code); printIndex(tindex, tsize); printf (" /");
  for (i = 0; i < tsize; i++) { printf(" "); printIndex(insTindex[i], tsize); }
  printf(" /");
  for (i = 0; i < tsize - 1; i++) { printf(" "); printIndex(delTindex[i], tsize); }
  printf("\n");
  for (i = 0; i < states; i++) { 
    if (currV[i] < INT_MIN/3) printf(" -inf/%2d", currTr[i]);
    else printf("%5d/%2d", currV[i], currTr[i]); 
    if ((i % 10) == 9) printf("\n"); 
  }
  printf("\n");
}
#endif /* __DEBUG */

inline int 
getFrame(int state, int tsize, int startoffset, int stopoffset)
     /* returns frame if <state> is coding, relies on startoffset,
        startlength, stopoffset and stoplength to be multiples of
        3, returns -1 if not coding */
{
  int f;
  int d = state - iStart - startoffset + 1;
  if (d < 0) return -1;
  if (state < (iStop + stopoffset)) return(d % 3);
  for (f = 1; f < tsize; f++) {
    if ((iInsAfter[(12 - f) % 3] + f) == state) return 0;
    if ((iInsAfter[(13 - f) % 3] + f) == state) return 1;
    if ((iInsAfter[(14 - f) % 3] + f) == state) return 2;
    if ((iDelAfter[(14 - f) % 3] + f - 1) == state) return 0;
    if ((iDelAfter[(15 - f) % 3] + f - 1) == state) return 1;
    if ((iDelAfter[(16 - f) % 3] + f - 1) == state) return 2;
  }
  return -1;
}

static int maxSize = 0;
static int *V  = NULL;
static int *tr = NULL;

int
Compute(const char *seq, int iPen, int dPen, int minScore, SV *rv, 
	int cdsIndex, int utrIndex, int startIndex, int stopIndex, int minLen)
{
  const char *p;
  int f, i, s, code;
  int iCurr, bPrev, bScore, maxScore, tmpPrev, tmpScore;
  signed char **cdsMat = mg[cdsIndex].m, *utrMat = mg[utrIndex].m[0];
  signed char **startPf = mg[startIndex].m, **stopPf = mg[stopIndex].m;

  /* initialize some more parameters */
  int tsize = mg[cdsIndex].order, tindex = 1, tableSize = 1;
  int startlen = mg[startIndex].frames, startoffset = mg[startIndex].offset;
  int stoplen = mg[stopIndex].frames, stopoffset = mg[stopIndex].offset;
  int *insTindex = (int *)malloc(sizeof(int) * tsize);
  int *delTindex = (int *)malloc(sizeof(int) * (tsize - 1));
  AV *av = NULL;

  /* allocate tables and compute the state indices  */
  int states = 2 + startlen + stoplen + 6*tsize, seqLen = strlen(seq); 
  int mSize = sizeof(int) * seqLen * states;
  int *currV, *prevV, *currTr;
  if (maxSize < mSize) { 
    maxSize = mSize; 
    V  = (int *)realloc(V,  maxSize); 
    tr = (int *)realloc(tr, maxSize);
  }
  currV  = V; currTr = tr;

  /* tindex will point to the position representing the last tsize chars on seq */
  for (i = 0; i < tsize; i++) tableSize *= 5;      /* size of score tables per frame */
  tindex = tableSize - 1;                          /* tindex now represents all N's */
  for (i = 0; i < tsize; i++) insTindex[i] = tindex; 
  for (i = 0; i < tsize - 1; i++) delTindex[i] = tindex;
  initIndices(tsize, startlen, stoplen);

#ifdef __DEBUG
  printInitStatus(cdsIndex, utrIndex, startIndex, stopIndex, states, seqLen, 
		  tsize, tableSize, tindex, insTindex, delTindex);
#endif

  /* fill in the Viterbi and traceback tables, initialize for first char on seq */
  code = GetCode(seq[0]); tindex = (5*tindex + code) % tableSize;
  currV[i5utr] = ts5uPen + utrMat[tindex];
  for (f = 0; f < startlen; f++) currV[iStart+f] = minScore + startPf[f][code];
  for (f = 0; f <  3; f++) currV[iCds+f] = tscPen + cdsMat[f][tindex];
  for (f = 0; f < stoplen; f++) currV[iStop+f] = minScore + stopPf[f][code];
  currV[i3utr] = ts3uPen + utrMat[tindex];
  for (s = i3utr + 1; s < states; s++) currV[s] = INT_MIN/2;
  for (s = 0; s < states; s++) currTr[s] = iBegin;

#ifdef __DEBUG
  printCurrentStatus(0, seq[0], code, tindex, tsize, 
		     insTindex, delTindex, states, V, tr);
#endif

  /* fill in the Viterbi and traceback tables, main part */
  for (p = seq + 1; *p; p++) {
    prevV = currV, currV += states, currTr += states;

    /* update index variables */
    code = GetCode(*p);
    for (i = tsize - 1; i > 0; i--) insTindex[i] = (5*insTindex[i-1] + code) % tableSize;
    for (i = tsize - 2; i > 0; i--) delTindex[i] = (5*delTindex[i-1] + code) % tableSize;
    insTindex[0] = tindex;
    delTindex[0] = (25*tindex + 20 + code) % tableSize;
    tindex = (5*tindex + code) % tableSize;

    /* consider current nucleotide in 5'UTR */
    /* transitions UTR->UTR and CDS->CDS are presumed zero */
    currV[i5utr]  = prevV[i5utr] + utrMat[tindex]; 
    currTr[i5utr] = i5utr;

    /* consider current nucleotide in start profile */
    currV[iStart]  = prevV[i5utr] + t5ucPen + startPf[0][code];
    currTr[iStart] = i5utr;
    iCurr = iStart; bPrev = iStart - 1;
    for (f = 1; f < startlen; f++) {
      iCurr++; bPrev++;
      currV[iCurr]  = prevV[bPrev] + startPf[f][code];
      currTr[iCurr] = bPrev;
    }
    
    /* consider current nucleotide in CDS */
    iCurr = iCds; bPrev = iCds - 1; bScore = prevV[bPrev];
    findMax(i5utr, prevV, minScore, &bPrev, &bScore);
    findMax0(iCurr + 2, prevV, &bPrev, &bScore);
    findMax0(iInsNext[0], prevV, &bPrev, &bScore);
    findMax0(iDelNext[0], prevV, &bPrev, &bScore);
    currV[iCurr] = bScore + cdsMat[0][tindex]; 
    currTr[iCurr] = bPrev;
    for (f = 1; f < 3; f++) {
      iCurr++; bPrev = iCurr - 1; bScore = prevV[bPrev];
      findMax0(iInsNext[f], prevV, &bPrev, &bScore);
      findMax0(iDelNext[f], prevV, &bPrev, &bScore);
      currV[iCurr] = bScore + cdsMat[f][tindex]; 
      currTr[iCurr] = bPrev;
    }

    /* consider current nucleotide in stop profile */
    bPrev = INT_MIN; bScore = INT_MIN;
    for (f = startoffset + 2; f < startlen; f += 3) 
      findMax(iStart + f, prevV, minScore, &bPrev, &bScore);
    for (f = 0; f < tsize; f++) 
      findMax(iInsAfter[(14 - f)%3] + f, prevV, minScore, &bPrev, &bScore);
    for (f = 0; f < tsize - 1; f++) 
      findMax(iDelAfter[(15 - f)%3] + f, prevV, minScore, &bPrev, &bScore);
    tmpPrev = bPrev; tmpScore = bScore;
    findMax(iCds + 2, prevV, tc3uPen, &bPrev, &bScore);
    currV[iStop]  = bScore + stopPf[0][code]; currTr[iStop] = bPrev;
    iCurr = iStop; bPrev = iStop - 1;
    for (f = 1; f < stoplen; f++) {
      iCurr++; bPrev++;
      currV[iCurr]  = prevV[bPrev] + stopPf[f][code];
      currTr[iCurr] = bPrev;
    }

    /* consider current nucleotide in 3' UTR */
    bPrev = tmpPrev; bScore = tmpScore;       /* jump stop profile */
    bPrev = INT_MIN; bScore = INT_MIN;        /* !!! possibly remov! */
    findMax(i3utr - 1, prevV, 0, &bPrev, &bScore);
    findMax0(i3utr, prevV, &bPrev, &bScore);
    findMax(iCds+2, prevV, minScore, &bPrev, &bScore);
    currV[i3utr]  = bScore + utrMat[tindex]; 
    currTr[i3utr] = bPrev;

    /* consider current nucleotide in CDS after insertion */
    for (f = 0; f < 3; f++) {
      iCurr = iInsAfter[f]; bPrev = iCds + f;
      currV[iCurr]  = prevV[bPrev] + iPen;
      currTr[iCurr] = bPrev; 
      iCurr = iInsAfter[f]; bPrev = iCurr - 1;
      for (i = 1; i < tsize; i++) {
	iCurr++; bPrev++;
	currV[iCurr] = prevV[bPrev] + cdsMat[(i+f)%3][insTindex[i]];
	currTr[iCurr] = bPrev;
      }
    }

    /* consider current nucleotide in CDS after deletion */
    for (f = 0; f < 3; f++) {
      iCurr = iDelAfter[f]; bPrev = iCds + f;
      currV[iCurr]  = prevV[bPrev] + dPen + cdsMat[(f+2)%3][delTindex[0]];
      currTr[iCurr] = bPrev; 
      iCurr = iDelAfter[f]; bPrev = iCurr - 1;
      for (i = 1; i < tsize - 1; i++) {
	iCurr++; bPrev++;
	currV[iCurr] = prevV[bPrev] + cdsMat[(i+f+2)%3][delTindex[i]];
	currTr[iCurr] = bPrev;
      }
    }

#ifdef __DEBUG
    printCurrentStatus(p-seq, *p, code, tindex, tsize, 
		       insTindex, delTindex, states, currV, currTr);
#endif
  }

  /* fill in the Viterbi and traceback tables, terminate and find best */
  bPrev = i5utr; bScore = currV[bPrev] + t5uePen; 
  for (f = 0; f < startlen; f++) findMax(iStart+f, currV, minScore, &bPrev, &bScore);
  for (f = 0; f < 3; f++) { findMax(iCds+f, currV, tcePen, &bPrev, &bScore);
    for (i = 0; i < tsize; i++)
      findMax(iInsAfter[f] + i, currV, tcePen, &bPrev, &bScore);
    for (i = 0; i < tsize - 1; i++)
      findMax(iDelAfter[f] + i, currV, tcePen, &bPrev, &bScore);
  }
  for (f = 0; f < stoplen; f++) findMax(iStop+f, currV, minScore, &bPrev, &bScore);
  findMax(i3utr, currV, t3uePen, &bPrev, &bScore);
#ifdef __DEBUG
  printf("finished to fill Viterbi matrix, best score %d in state %d\n",bScore,bPrev);
#endif

  /* traceback and generate coding sequences starting from bPrev (confidence bScore) */
  maxScore = INT_MIN;
  if ( SvROK( rv ) && (SvTYPE(SvRV(rv)) == SVt_PVAV) ) av = (AV*)SvRV(rv);
  iCurr = bPrev; p--;
  while(iCurr != iBegin) {
    char *res = (char *)malloc(sizeof(char) * 2 * seqLen);
    AV *resArray = newAV(); 
    SV *sv; 
    int iOld = iCurr, rStart, rStop;
    char  *r, *q;

    /* skip non coding */
    while(((iBegin < iCurr)
	   && (iCurr < iStart + startoffset - 1))
	  || (((iStop + stopoffset - 1) < iCurr)
	      && (iCurr < iInsAfter[0]))) { 
#ifdef __DEBUG
      printf("trace back non-coding: state %d position %4d(%c)\n", iCurr, (p-seq), *p);
#endif
      iCurr = currTr[iCurr]; 
      currTr -= states;
      p--; 
    }

    /* treat coding */
    if (iCurr != iBegin) {
      int rScore = V[(p-seq)*states + iCurr];
      r = res;
      rStop = (p-seq);
      if (getFrame(iCurr, tsize, startoffset, stopoffset) == 0) {*r++ = 'X';*r++ = 'X';}
      if (getFrame(iCurr, tsize, startoffset, stopoffset) == 1) *r++ = 'X';
      while(((iStart + startoffset - 1 <= iCurr)
	     && (iCurr <= iStop + stopoffset - 1))
	    || (iInsAfter[0] <= iCurr)) {
	int done = 0;
	for (f = 0; f < 3; f++) {
	  if (iCurr == iInsAfter[f]) { *r++ = tolower(*p); done = 1; }
	  if (iCurr == iDelAfter[f]) { *r++ = toupper(*p); *r++ = 'X'; done = 1; }
	}
	if (done == 0) *r++ = toupper(*p);
	/* remove stop-profile penalty from coding score */
	if ((iCurr == iCds + 2) && (iOld == iStop)) rScore -= tc3uPen;
#ifdef __DEBUG
      printf("trace back     coding: state %2d(%2d) position %4d(%c)\n", 
	     iCurr, getFrame(iCurr, tsize, startoffset, stopoffset), p-seq, *(r-1));
#endif
	iOld = iCurr;
	iCurr = currTr[iCurr];
	currTr -= states;
	p--;
      }
      rStart = p-seq + 1;
      if (p >= seq) rScore -= V[(p-seq)*states + iCurr];
      if (rScore > maxScore) maxScore = rScore;
      if (getFrame(iOld, tsize, startoffset, stopoffset) == 1) *r++ = 'X';
      if (getFrame(iOld, tsize, startoffset, stopoffset) == 2) { *r++='X'; *r++='X'; }
      *r-- = 0;

      /* reverse the array and add to the result-array */
      q = res;
      while (q < r) { char c = *r; *r-- = *q; *q++ = c; } 
#ifdef __DEBUG
      printf("found coding %s, add to results, state %d \n", res, iCurr);
#endif
      if ( SvROK( rv ) && (SvTYPE(SvRV(rv)) == SVt_PVAV) && (rStop - rStart >= minLen)) {
	sv = newSViv(rScore); av_push(resArray, sv);
	sv = newSViv(rStart); av_push(resArray, sv);
	sv = newSViv(rStop);  av_push(resArray, sv);
	sv = newSVpv(res, 0); av_push(resArray, sv);
	/* Create a reference to the new array, and push it on the results array.  */ 
	sv = newRV_inc((SV*) resArray); av_push(av, sv);
      }
    }
    free(res);
  }

  /* clean up */
  free(insTindex);
  free(delTindex);

  return maxScore;
}
