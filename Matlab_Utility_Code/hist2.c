/*

  HIST2.C	fast two-dimensional histogramming code

  [C, BX, BY] = HIST2 (DATA, BINX, BINY).

  Calculate 2D histogram.  C is the counts in each bin (matrix).  BX is
  the vector of X centers, BY is the matrix of Y centers.

  DATA is ndata by 2.

  BINX, BINY are vectors of bin boundaries (NOTE: this is different to
  the mathworks HIST which takes bin centers) If BINi is a scalar it
  is taken to be a bin count.

  If both BINX and BINY are specified as a count, or as a regularly
  spaced vector, a rapid binning algorithm is used that depends on
  this and is linear in ndata.  If BINi is irregularly spaced the
  algorithm used takes time proportional to ndata*nbins.

*/

#include <stdio.h>
#include <math.h>
#include "mex.h"

/* Arguments */

#define	DATA_IN	   prhs[0]
#define	BINX_IN	   prhs[1]
#define	BINY_IN	   prhs[2]
#define	C_OUT	   plhs[0]
#define	BX_OUT	   plhs[1]
#define	BY_OUT	   plhs[2]

/* Auxilliary prototypes */
static void findext (double *, int, double *, double *);


/******************************************************************************
  WORKHORSE FUNCTIONS
  */

/*
 * regular bin size histogramming routine (faster)
 */

static void reghist2d(double	 data[], 
		      int	 ndata, 
		      double	 minx, 
		      double	 sizex, 
		      int	 nbinx,
		      double	 miny, 
		      double	 sizey, 
		      int	 nbiny,
		      double	 cnts[],
		      double	 ctrx[],
		      double	 ctry[])
{
  int i, bx, by;
  double maxx = minx + sizex * nbinx;
  double maxy = miny + sizey * nbiny;

#define datax(i) data[i]
#define datay(i) data[ndata + i]

  for (i = 0; i < nbinx*nbiny; i++) {
    cnts[i] = 0;
  }

  for (i = 0; i < nbinx; i++) {
    ctrx[i] = minx + i*sizex + sizex/2;
  }

  for (i = 0; i < nbiny; i++) {
    ctry[i] = miny + i*sizey + sizey/2;
  }

  for (i = 0; i < ndata; i++) {
    if (datax(i) < minx || datax(i) >= maxx || 
	datay(i) < miny || datay(i) >= maxy)
      continue;
    bx = (int)((datax(i)-minx)/sizex);
    by = (int)((datay(i)-miny)/sizey);
    if (bx < nbinx && by < nbiny)
      cnts[bx * nbiny + by] ++;
  }

#undef datax
#undef datay
}


static void binhist(double	 data[], 
		    int		 ndata, 
		    double	 bins[],
		    int		 nbins,
		    double	 cnts[],
		    double	 ctrs[])
{
  int i, j;

  for (i = 0; i < nbins; i++) {
    cnts[i] = 0;
    ctrs[i] = (bins[i] + bins[i+1])/2;
  }

  for (i = 0; i < ndata; i++)
    for (j = 0; j < nbins; j++)
      if (data[i] >= bins[j] && data[i] < bins[j+1]) 
	cnts[j] ++;
}


/******************************************************************************
  INTERFACE FUNCTION
  */

void mexFunction(int		nlhs,
		 Matrix	*	plhs[],
		 int		nrhs,
		 Matrix	*	prhs[])
{
  double	*data = NULL, *binx = NULL, *biny = NULL;
  double	*cnts = NULL, *ctrx = NULL, *ctry = NULL;
  int		ndata, nbinx, nbiny;
  double 	minx, maxx, sizex = -1, miny, maxy, sizey = -1;
  double	*t,*y;
  unsigned int	i,m,n;

  /* Check numbers of arguments */
  if (nrhs == 0) {
    mexErrMsgTxt("HIST2D: no data to histogram");
  } else if (nrhs > 3) {
    mexErrMsgTxt("HIST2D: two many arguments.");
  }
  if (nlhs < 3) {
    mexErrMsgTxt("HIST2D: must be called with three output arguments");
  }

  /* Get data */
  m = mxGetM(DATA_IN);
  n = mxGetN(DATA_IN);
  ndata = m;
  if (!mxIsNumeric(DATA_IN) || mxIsComplex(DATA_IN) || 
      !mxIsFull(DATA_IN)  || !mxIsDouble(DATA_IN)) {
    mexErrMsgTxt("HIST2D: data must be a full real valued matrix.");
  }
  if (n != 2)
    mexErrMsgTxt("HIST2D: data must have two columns");

  data = mxGetPr(DATA_IN);


  /* Get bin specification */

  /* x bins */
  m = mxGetM(BINX_IN);
  n = mxGetN(BINX_IN);
  if (!mxIsNumeric(BINX_IN) || mxIsComplex(BINX_IN) ||
      !mxIsFull(BINX_IN) || !mxIsDouble(BINX_IN) ||
      (m != 1 && n != 1)) {
    mexErrMsgTxt("HIST2D: bin spec must be a real scalar or vector");
  }
  if (m == 1 && n == 1) {	/* number of bins */
    nbinx = (int)*(double *)mxGetPr(BINX_IN);
    findext (data, ndata, &minx, &maxx);
    sizex = (maxx-minx)/nbinx;
  } else {			/* vector of bin boundaries */
    nbinx = n*m - 1;
    binx  = mxGetPr (BINX_IN);

    /* check to see if spacing is regular -- if so use fast algorithm */
    sizex = binx[1] - binx[0];
    for (i = 1; i < nbinx; i++)
      if (fabs(binx[i+1] - binx[i] - sizex) > 1e-3*sizex) {
	sizex = 0;
	break;
      }
    if (sizex) {
      minx = binx[0];
    }
  }


  /* y bins */
  m = mxGetM(BINY_IN);
  n = mxGetN(BINY_IN);
  if (!mxIsNumeric(BINY_IN) || mxIsComplex(BINY_IN) ||
      !mxIsFull(BINY_IN) || !mxIsDouble(BINY_IN) ||
      (m != 1 && n != 1)) {
    mexErrMsgTxt("HIST2D: bin spec must be a real scalar or vector");
  }
  if (m == 1 && n == 1) {	/* number of biny */
    nbiny = (int)*(double *)mxGetPr(BINY_IN);
    findext (data+ndata, ndata, &miny, &maxy);
    sizey = (maxy-miny)/nbiny;
  } else {			/* vector of bin boundaries */
    nbiny = n*m - 1;
    biny  = mxGetPr (BINY_IN);

    /* check to see if spacing is regular -- if so use fast algorithm */
    sizey = biny[1] - biny[0];
    for (i = 1; i < nbiny; i++)
      if (fabs(biny[i+1] - biny[i] - sizey) > 1e-3*sizey) {
	sizey = 0;
	break;
      }
    if (sizey) {
      miny = biny[0];
    }
  }


    
  /* Create output matrices */
  C_OUT = mxCreateFull(nbiny, nbinx, REAL);
  cnts = mxGetPr(C_OUT);

  BX_OUT = mxCreateFull(1, nbinx, REAL);
  ctrx = mxGetPr(BX_OUT);
  
  BY_OUT = mxCreateFull(1, nbiny, REAL);
  ctry = mxGetPr(BY_OUT);

  /* Do the actual computations in a subroutine */

  if (sizex && sizey) {
    reghist2d(data, ndata, 
	      minx, sizex, nbinx, 
	      miny, sizey, nbiny,
	      cnts, ctrx, ctry);
  } else {
    mexErrMsgTxt ("HIST2D: irregular binsize histogramming not implemented");
    /*    binhist2d(data, ndata, 
	      binx, nbinx, 
	      biny, nbiny, 
	      cnts, ctrx, ctry);
	      */
  }

  return;
}


/******************************************************************************
  AUXILLIARY FUNCTIONS
  */


static void findext(double	 data[], 
		    int		 ndata, 
		    double	*min, 
		    double      *max)
{
  int i;

  *min = *max = data[0];
  for (i = 1; i < ndata; i++) {
    if (data[i] < *min) *min = data[i];
    if (data[i] > *max) *max = data[i];
  }
}
