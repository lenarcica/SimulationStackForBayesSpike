#ifndef HEAPSORTH
  #define HEAPSORTH 1
  #include "HeapSort.h"
#endif
extern "C" {
SEXP HeapSortShell(SEXP SDoubles, SEXP OutII) {
  int *Out = HeapSortAll(length(SDoubles), REAL(SDoubles));
  int ii;
  if (isReal(OutII)) {
   for (ii = 0; ii < length(OutII); ii++) {
     REAL(OutII)[ii] = (double) Out[ii];
   }
  } else if (isInteger(OutII)) {
   for (ii = 0; ii < length(OutII); ii++) {
     INTEGER(OutII)[ii] = (int) Out[ii];
   }  
  }
  Free(Out);
  return(OutII);  
}
}

int *HeapSortAll(int N, double *Doubles) {
int *OutIntii = (int*)Calloc(N, int);
unsigned long i=0, ir, j, l;
double rra = 0.0;

int ii = 0;  int rraii = 0;
for (ii = 0; ii < N; ii++) { OutIntii[ii] = ii; }

if (N < 2) return( OutIntii );
l = (N >> 1) +1;  ir = N;
for (;;) {
  if (l > 1) {
    l--;
    rra = Doubles[ OutIntii[l-1] ];  rraii = OutIntii[l-1];
  } else {
    rra =  Doubles[ OutIntii[ir-1] ];  rraii = OutIntii[ir-1];
    OutIntii[ir-1] = OutIntii[0];
    if (--ir == 1) {
      OutIntii[0] = rraii;     // Least Competent worker of all?
      break;
    }
  }
  i=l; j = l + l;
  while( j <= ir) {
    if (j < ir && Doubles[OutIntii[j-1]] < Doubles[OutIntii[j]]) j++;
    if (rra < Doubles[OutIntii[j-1]]) {
      OutIntii[i-1] = OutIntii[j-1];
      i = j;
      j <<= 1;
    }  else {
      break;
    } 
  }
  OutIntii[i-1] = rraii;
}
return(OutIntii);
}

int ReverseString(int N, int *ReverseScores) {
  int HalfN = floor(N/2.0);
  int ii;
  for (ii = 0; ii < HalfN; ii++) {
    if (N-ii-1 >= N) {
      Rprintf("ReverseString: Absoulte failure N-ii-1 = ", N-ii-1);
      R_FlushConsole();
    }
    ReverseScores[ii] += ReverseScores[N-ii-1];
    ReverseScores[N-ii-1] = ReverseScores[ii] - ReverseScores[N-ii-1];
    ReverseScores[ii] = ReverseScores[ii] - ReverseScores[N-ii-1];  
  }
  return(1);
}