/*********************************************************************
 * Demo.cpp
 *
 * This file shows the basics of setting up a mex file to work with
 * Matlab.  This example shows how to use 2D matricies.  This may
 * 
 * Keep in mind:
 * <> Use 0-based indexing as always in C or C++
 * <> Indexing is column-based as in Matlab (not row-based as in C)
 * <> Use linear indexing.  [x*dimy+y] instead of [x][y]
 *
 * For more information, see my site: www.shawnlankton.com
 * by: Shawn Lankton
 *
 ********************************************************************/
#include <matrix.h>
#include "GCoptimization.h"
#include <mex.h>   
#include <stdlib.h>
/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//declare variables
    mxArray *image_in_m, *mask_in_m, *c_out_m;
    const mwSize *dims;
    double *image, *mask, *c, sigma;
    int MXS, MYS, numdims;
    int x,y;

//associate inputs
    image_in_m = mxDuplicateArray(prhs[0]);
    mask_in_m = mxDuplicateArray(prhs[1]);
    //d_in_m = mxDuplicateArray(prhs[2]); 
//figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    MYS = (int)dims[0]; MXS = (int)dims[1];

//associate outputs
    c_out_m = plhs[0] = mxCreateDoubleMatrix(MYS,MXS,mxREAL);
   
//associate pointers
    image = mxGetPr(image_in_m);
    mask = mxGetPr(mask_in_m);
    sigma = mxGetScalar(prhs[2]);
    c = mxGetPr(c_out_m);

//do something
    GCoptimizationGeneralGraph gc(MXS*MYS);
   /*
    for(y=0;y<MYS;y++)
    {
        for(x=0;x<MXS;y++)
        {
           
            c[i*dimy+j] = a[i*dimy+j]+d; //adds 5 to every element in a
            mexPrintf("element[%d][%d] = %f\n",j,i,c[i*dimy+j]);
        }
    }
*/
    return;
}