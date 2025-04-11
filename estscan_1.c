/* $Id: estscan_1.c,v 1.1.1.1 2006/12/18 15:44:47 c4chris Exp $
 *
 * Christian Iseli, LICR ITO, Christian.Iseli@licr.org
 *
 * Copyright (c) 1999 Swiss Institute of Bioinformatics. All rights reserved.
 */
#include "estscan_1.h"
#include <stdio.h>
#include <limits.h>
#ifndef __GNUC__
#include <alloca.h>
#endif

/* At this point, we do not allow the deletion of unused matrices...  */

typedef struct _matrixGroup {
  signed char **m;
  int order;
  int frames;
  int offset;
} matrixGroup_t;

#define ALLOC_CHUNK 10

static matrixGroup_t *mg;
static int groupMax;
static int nextGroup;

static int
intCompare(const void *a, const void *b)
{
  return (*(int *)a < *(int *)b) ? -1 : ((*(int *)a > *(int *)b) ? 1 : 0);
}

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

int
CreateMatrix(int order, int frames, int offset, SV *rv, int min, int method)
{
  AV *av;
  int avlen, i, frame;
  int sSize = 1;
  int *step = (int *)alloca(sizeof(int) * order);
  int *sStep = (int *)alloca(sizeof(int) * order);
  signed char **fScore = (signed char **)alloca(sizeof(signed char *) * frames);

  /*printf("%d. call to CreateMatrix: %d, %d, %d, %d, %d\n", nextGroup, order, frames, offset, min, method); !!!*/

  if ( SvROK( rv ) && (SvTYPE(SvRV(rv)) == SVt_PVAV) )
    av = (AV*)SvRV(rv);
  else {
    warn("CreateMatrix: 2nd parameter was not an array ref");
    return -1;
  }
  if (order < 1) {
    warn("CreateMatrix: order should be >=1 (%d)", order);
    return -1;
  }
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
  for (i = 0; i < order; i++)
    sSize *= 5;
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
		die("Bad method (%d) to compute N score value.", method);
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
  mg[nextGroup].order = order;
  mg[nextGroup].frames = frames;
  mg[nextGroup].offset = offset;
  return nextGroup++;
}

static int
GetCode(char c)
{
  c &= 0x1f; /* Get lower bits, get rid of upper/lower info.  */
  switch (c) {
  case 1: /* This is A.  */
    return 0;
  case 3: /* This is C.  */
    return 1;
  case 7: /* This is G.  */
    return 2;
  case 20: /* This is T.  */
    return 3;
  default: /* Everything else.  */
    return 4;
  }
}

#if 0 /* CI */
static void
printIndex(int index)
{
  char s[7];
  int i;
  s[6] = 0;
  for (i = 5; i >= 0; i--) {
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
  printf("%s", s);
}
#endif /* 0 CI */

static void
reconstruct(char *seq, char *trace[3], int tPos, char *res, int frame, int *start)
{
  char *s = seq + tPos;
  char *r = res;
  int done = 0;
  while (s >= seq && !done) {
    switch (trace[frame][tPos]) {
    case 'e':
      *r++ = toupper(*s);
      s -= 1;
      tPos -= 1;
      frame = (frame == 0) ? 2 : frame - 1;
      break;
    case 'd':
      *r++ = toupper(*s); s -= 1;
      *r++ = toupper(*s); s -= 1;
      *r++ = toupper(*s); s -= 1;
      *r++ = toupper(*s); s -= 1;
      *r++ = toupper(*s); s -= 1;
      *r++ = tolower(*s); s -= 1;
      frame = (frame == 2) ? 0 : frame + 1;
      tPos -= 6;
      break;
    case 'i':
      *r++ = toupper(*s); s -= 1;
      *r++ = toupper(*s); s -= 1;
      *r++ = toupper(*s); s -= 1;
      *r++ = toupper(*s); s -= 1;
      *r++ = toupper(*s); s -= 1;
      *r++ = toupper(*s); s -= 1;
      *r++ = 'X';
      frame = (frame == 0) ? 2 : frame - 1;
      tPos -= 6;
      break;
    case 'b':
      if (*s) *r++ = toupper(*s); s -= 1;
      if (*s) *r++ = toupper(*s); s -= 1;
      if (*s) *r++ = toupper(*s); s -= 1;
      if (*s) *r++ = toupper(*s); s -= 1;
      if (*s) *r++ = toupper(*s); s -= 1;
      frame = (frame == 2) ? 0 : frame + 1;
      done = 1;
      break;
    case 's':
      done = 1;
      break;
    default:
      die("Bad letter %c in trace.", *s);
    }
  }
  /* The following seems to manage to determine the correct reading
     frame.  */
  switch (frame) {
  case 1:
    if (s >= seq) {
      *r++ = toupper(*s); s -= 1;
    } else 
      *r++ = 'X';
    /* Fall through.  */
  case 0:
    if (s >= seq) {
      *r++ = toupper(*s); s -= 1;
    } else 
      *r++ = 'X';
  default:
    /* Do nothing.  */;
  }
  /* Determine starting position.  */
  *start = s - seq + 1;
  *r-- = 0;
  while (res < r) {
    char c = *r;
    *r-- = *res;
    *res++ = c;
  }
}

static int
evaluateProfileAt(char *seq, int mindex, int index) 
  /* Evaluates the profile pointed to by mindex at 
     position index and returns the corresponding score,
     uncomplete evaluation returns -100 */
{
  int i;
  int score = 0;
  char *ptr;
  if (index < mg[mindex].offset) { return -100; }
  ptr = seq + index  - mg[mindex].offset; 
  if (strlen(ptr) < mg[mindex].frames) { return -100; }
  for (i = 0; i < mg[mindex].frames; i++) {
    /*    printf("%c", *ptr); */
    score += mg[mindex].m[i][GetCode(*ptr++)];
  }
  /*  printf(" yields score %d\n", score);*/
  return score;
}

static void
finetune_stop(char *seq, int maxPos, int stopIndex, int stopDist, char *res, int start, int *end)
     /* within the distance stopDist from *end, searches
	for end codons aligned with the frame given by start.
	(stopIndex is currently not used, instead codons are
	searched for explicitely). While searching deletions
        (small characters) are skipped. The first stop/codon
	encountered is interpreted as end of the coding 
	sequece. res and *end are modified accordingly. 
	maxPos indicates where the non-coding part of seq
	starts.  */
{
  int currFrame;
  char *ptr, *from, *to, c[3];

  /*printf("finetune seq %s,\nmaxPos %d, stopDist %d, start %d, end %d\nres %s\n", 
	 seq, maxPos, stopDist, start, *end, res);*/

  /* scan into window determined by stopDist */
  ptr = res;
  from = res + strlen(res) - stopDist;
  currFrame = 0;
  while ((ptr < from) || (currFrame != 0)) { 
    /*printf("skip %c in frame %d\n", *ptr, currFrame);*/
  if (*ptr < 'Z') currFrame = (currFrame + 1) % 3;
    ptr++; 
  }

  /* search stop codons within res */
  while (*ptr > 'Z') ptr++; 
  c[0] = *ptr++;
  currFrame = 1;
  while (*ptr) {
    while (*ptr > 'Z') ptr++; 
    c[currFrame] = *ptr;
    currFrame = (currFrame + 1) % 3;
    ptr++; 
    if (currFrame == 0) {
      /*printf("In res, is %c%c%c a stop codon?\n", c[0], c[1], c[2]);*/
      if (((c[0] == 'T') && (c[1] == 'A') && (c[2] == 'A')) ||
	  ((c[0] == 'T') && (c[1] == 'G') && (c[2] == 'A')) ||
	  ((c[0] == 'T') && (c[1] == 'A') && (c[2] == 'G'))) { 
	ptr -= 3;
	/*printf("fine tuned stop from %d, ptr points to %c ", *end, *ptr);*/
	*ptr++ = 0; *end -= 1;
	while (*ptr++) *end -= 1;
	/*printf("to %d due to %c%c%c\n", *end, c[0], c[1], c[2]);*/
	return;
      }
    }
  }

  /* search stop codons behind res */
  ptr = seq + maxPos + 1;
  if (strlen(ptr) < stopDist) { to = ptr + strlen(ptr); }
  else { to = ptr + stopDist; }
  while (ptr < to - 2) {
    c[0] = *ptr++;
    c[1] = *ptr++;
    c[2] = *ptr++;
    /*printf("After res, is %c%c%c a stop codon\n", c[0], c[1], c[2]);*/
    if (((c[0] == 'T') && (c[1] == 'A') && (c[2] == 'A')) ||
	((c[0] == 'T') && (c[1] == 'G') && (c[2] == 'A')) ||
	((c[0] == 'T') && (c[1] == 'A') && (c[2] == 'G'))) { 
      ptr -= 3;
      /*printf("fine tuned stop from %d ", *end);*/
      res += strlen(res);
      seq += maxPos + 1;
      while (seq < ptr) { *res++ = *seq++; *end += 1; }
      *res = 0;       
      /*printf("to %d due to %c%c%c\n", *end,  c[0], c[1], c[2]);*/
      return;
    }
  }
}

typedef struct _sInfo_t {
  int index;
  int score;
} sInfo_t;

#define INSERTION 0
#define DELETION 1

int
Compute(char *seq, int iPen, int dPen, int cutOff, int dropMax, SV *rv, int gIndex, 
	int startIndex, int stopIndex, int stopDist, int ncIndex, char *graph)
{
  int tupleSize = mg[gIndex].order;
  AV *av = NULL;
  int sLen = strlen(seq);
  char *ptr = seq;
  int startPos = 0;
  int bigMax = 0; /* The overall maximum.  */
  int maxOnly = 0;
  int max = (startIndex == -1) ? 0 : INT_MIN;
  int maxPos = 0;
  int maxFrame = 0;
  char *res = malloc(sizeof(char) * (sLen * 2 + stopDist));
  sInfo_t **wing[2];
  sInfo_t mLine[3];
  int E, f, i;
  signed char *fScore[3];
  unsigned int tableMod = 1;
  FILE *codingfh = NULL;
  FILE *startfh = NULL;
  char *trace[3]; /* Easier to pass to reconstruct.  */
  trace[0] = malloc(sizeof(char) * (sLen + 1));
  trace[1] = malloc(sizeof(char) * (sLen + 1));
  trace[2] = malloc(sizeof(char) * (sLen + 1));

  if (sLen < tupleSize) return 0;
  if ( SvROK( rv ) && (SvTYPE(SvRV(rv)) == SVt_PVAV) ) av = (AV*)SvRV(rv);
  else maxOnly = 1;
  if (gIndex < 0 || gIndex >= nextGroup) die("Bad group index (%d)", gIndex);
  wing[0] = alloca(sizeof(sInfo_t *) * tupleSize);
  wing[1] = alloca(sizeof(sInfo_t *) * tupleSize);
  for (i = 0; i < tupleSize; i++) {
    tableMod *= 5;
    wing[0][i] = alloca(sizeof(sInfo_t) * 3);
    wing[1][i] = alloca(sizeof(sInfo_t) * 3);
  }

  if (graph[0]) { 
    char *filename = alloca(sizeof(char) * strlen(graph) + 10);
    strcpy(filename, graph);
    strcat(filename, ".coding");
    codingfh = fopen(filename, "w"); 
    if (startIndex > -1) {
      strcpy(filename, graph);
      strcat(filename, ".start");
      startfh = fopen(filename, "w"); 
    }
  }

  for (i = 0; i < 3; i++) fScore[i] = mg[gIndex].m[i];
  while (*ptr != 0) {
    int curPos = tupleSize - 1;
    /* Set things up.  */

    max = (startIndex == -1) ? 0 : INT_MIN;
    maxPos = 0;
    maxFrame = 0;
    for (f = 0; f < 3; f++) {
      if (maxOnly == 0) memset(trace[f], 'e', tupleSize);
      for (i = 0; i < 2; i++) {
	for (E = 0; E < tupleSize; E++) {
	  wing[i][E][f].index = 0;
	  wing[i][E][f].score = -1000;
	}
      }
    }
    /* Compute first indices and scores.  */  
    mLine[0].index = 0;
    for (i = 0; i < tupleSize; i++) mLine[0].index = mLine[0].index * 5 + GetCode(*ptr++);
    mLine[0].score = fScore[0][mLine[0].index];
    for (f = 1; f < 3; f++) {
	mLine[f].index = mLine[0].index;
	mLine[f].score = fScore[f][mLine[0].index];
    }
    /* Determine the best main score.  */
    for (f = 0; f < 3; f++) {
      if ((mLine[f].score < 0) && (startIndex == -1)){
	mLine[f].score = 0;
	trace[f][tupleSize - 1] = 'b';
      }
      if (mLine[f].score > max) {
	max = mLine[f].score;
	if (max > bigMax) bigMax = max;
	maxPos = tupleSize - 1;
	maxFrame = f;
      }
    }
    if (codingfh) { fprintf(codingfh, "%d %d\n", curPos + startPos, max); }
    while (*ptr != 0) {
      int index;
      int curC = GetCode(*ptr++);
      int curScore[3];
      curPos += 1;
      /* Compute the next scores.  The index is only updated in frame 0.  */
      for (i = 0; i < 2; i++) {
	for (E = 0; E < tupleSize - 1; E++) {
	  /* Simple shifting for the wings.  */
	  index = (wing[i][E + 1][0].index * 5) % tableMod + curC;
	  wing[i][E][0].index = index;
	  for (f = 0; f < 3; f++) {
	    int prevF = (f == 0) ? 2 : f - 1;
	    wing[i][E][f].score = wing[i][E + 1][prevF].score + fScore[f][index];
	  }
	}
      }
      /* Now for the interesting stuff.  First, the delete branch.  */
      /* We keep the current index.  */
      wing[DELETION][tupleSize-1][0].index = mLine[0].index;
      for (f = 0; f < 3; f++) {
	/* No frame shift, no score to add, except for the deletion penalty.  */
	wing[DELETION][tupleSize-1][f].score = mLine[f].score + dPen;
      }
      /* Now, the main branch.  First, keep a copy of the score.  */
      for (f = 0; f < 3; f++) curScore[f] = mLine[f].score;
      index = (mLine[0].index * 5) % tableMod + curC;
      mLine[0].index = index;
      for (f = 0; f < 3; f++) {
	int prevF = (f == 0) ? 2 : f - 1;
	mLine[f].score = curScore[prevF] + fScore[f][index];
      }
      /* And finally, the insertion branch, which can reuse the new main branch.  */
      index = (mLine[0].index * 5) % tableMod + 4;
      wing[INSERTION][tupleSize-1][0].index = index;
      for (f = 0; f < 3; f++) {
	int prevF = (f == 0) ? 2 : f - 1;
	wing[INSERTION][tupleSize-1][f].score = mLine[prevF].score + fScore[f][index] + iPen;
      }
      /* Determine the new best main score  */
      for (f = 0; f < 3; f++) {
	trace[f][curPos] = 'e';
	if (wing[DELETION][0][f].score > mLine[f].score) {
	  mLine[f].score = wing[DELETION][0][f].score;
	  /* index should be fine.  */
	  trace[f][curPos] = 'd';
	}
	if (wing[INSERTION][0][f].score > mLine[f].score) {
	  mLine[f].score = wing[INSERTION][0][f].score;
	  /* index should be fine.  */
	  trace[f][curPos] = 'i';
	}
	if ((mLine[f].score < 0) && (startIndex == -1)) {
	  mLine[f].score = 0;
	  trace[f][curPos] = 'b';
	}
	if (mLine[f].score > max) {
	  max = mLine[f].score;
	  if (max > bigMax) bigMax = max;
	  maxPos = curPos;
	  maxFrame = f;
	}
      }
      /* evaluate start model if required */
      if (startIndex > -1) {
	int startScore = evaluateProfileAt(seq, startIndex, curPos);
	if (startfh) fprintf(startfh, "%d %d\n", curPos, startScore - 100);
	if (mLine[0].score < startScore) {
	  /*printf("move start to %d with score %d, instead of %d, start %c%c%c\n", curPos, 
		 startScore, mLine[0].score, seq[curPos], seq[curPos+1], seq[curPos+2]); */
	  mLine[0].score = startScore;
	  trace[0][curPos] = 's';
	}
      }
      /* See if we have reached a maximum.  */
      if (max >= cutOff) {
	int notYet = 0;
	for (f = 0; f < 3; f++) 
	  if (mLine[f].score > 0 && max - mLine[f].score < dropMax)  notYet = 1;
	if (notYet == 0) {
	  if (maxOnly == 0) {
	    int seqStart, seqEnd;
	    AV *resArray = newAV();
	    SV *sv;
	    /* We have found a max, and are ready to restart.  */
	    while (maxFrame != 2 && maxPos + startPos < sLen) {
	      maxFrame += 1;
	      trace[maxFrame][++maxPos] = 'e';
	    }
	    reconstruct(seq + startPos, trace, maxPos, res, maxFrame, &seqStart);
	    seqStart += startPos;
	    seqEnd = maxPos + startPos;
	    if (stopIndex > -1) 
	      finetune_stop(seq, maxPos, stopIndex, stopDist, res, seqStart, &seqEnd);
	    /* Put the resulting sequence in the results array.  */
	    /* We create an array containing max, position, sequence.  */
	    sv = newSViv(max);
	    av_push(resArray, sv);
	    sv = newSViv(seqStart);
	    av_push(resArray, sv);
	    sv = newSViv(seqEnd);
	    av_push(resArray, sv);
	    sv = newSVpv(res, 0);
	    av_push(resArray, sv);
	    /* Create a reference to the new array, and push it on the results array.  */
	    sv = newRV_inc((SV*) resArray);
	    av_push(av, sv);
	  }
	  /* Restart the big loop.  */
	  startPos += maxPos + 1;
	  ptr = seq + startPos;
	  break;
	}
      }

      if (codingfh) { 
	int bestFrame = 0;
	if (mLine[0].score > bestFrame) { bestFrame = mLine[0].score; }
	if (mLine[1].score > bestFrame) { bestFrame = mLine[1].score; }
	if (mLine[2].score > bestFrame) { bestFrame = mLine[2].score; }
	fprintf(codingfh, "%d %d\n", curPos + startPos, bestFrame); 
      }
      /* CI printf("Pos: %d, scores: %d, %d, %d; max: %d\n", curPos + startPos,
	       mLine[0].score, mLine[1].score, mLine[2].score, max); */
    }
  }

  /* return results */
  if (max >= cutOff && maxOnly == 0) {
    int seqStart, seqEnd;
    AV *resArray = newAV();
    SV *sv;
    while (maxFrame != 2 && maxPos + startPos < sLen) {
      maxFrame += 1;
      trace[maxFrame][++maxPos] = 'e';
    }
    reconstruct(seq + startPos, trace, maxPos, res, maxFrame, &seqStart);
    seqStart += startPos;
    seqEnd = maxPos + startPos;
    if (stopIndex > -1) 
      finetune_stop(seq, maxPos, stopIndex, stopDist, res, seqStart, &seqEnd);
    /* Put the resulting sequence in the results array.  */
    /* We create an array containing max, position, sequence.  */
    sv = newSViv(max);
    av_push(resArray, sv);
    sv = newSViv(seqStart);
    av_push(resArray, sv);
    sv = newSViv(seqEnd);
    av_push(resArray, sv);
    sv = newSVpv(res, 0);
    av_push(resArray, sv);
    /* Create a reference to the new array, and push it on the results array.  */
    sv = newRV_inc((SV*) resArray);
    av_push(av, sv);
  }

  /* clean up and return */
  if (codingfh) { fprintf(codingfh, "\n\n"); fclose(codingfh); }
  if (startfh) { fprintf(startfh, "\n\n"); fclose(startfh); }
  free(trace[0]);
  free(trace[1]);
  free(trace[2]);
  free(res);
  return bigMax;
}
