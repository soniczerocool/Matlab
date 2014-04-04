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
#include <cmath>
#include <math.h>
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
    mxArray *image_in_m,  *c_out_m;
    const mwSize *dims;
    double *image, *c, beta;
    int MXS, MYS, numdims;
    int x,y;
    
    /* Check for proper number of input and output arguments */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt( "MATLAB: markov_network_multilable:invalidNumInputs",
                "Wrong number of arguments\n\n\t Data usage : labelfield = demo2(image,mask,sigma).");
    }
    
    /* Check data type of input argument  */
    if (!(mxIsDouble(prhs[0])) || !(mxIsDouble(prhs[1]))){
        mexErrMsgIdAndTxt( "MATLAB: markov_network_multilable:inputNotDouble",
                " Input arguments must be of type double.");
    }
    
    
//associate inputs
    image_in_m = mxDuplicateArray(prhs[0]);
 
    //d_in_m = mxDuplicateArray(prhs[2]); 
//figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    MYS = (int)dims[0]; MXS = (int)dims[1];

//associate outputs
    c_out_m = plhs[0] = mxCreateDoubleMatrix(MYS,MXS,mxREAL);
    
   
//associate pointers and scalars
    image = mxGetPr(image_in_m);
   
    beta = mxGetScalar(prhs[1]);
    c = mxGetPr(c_out_m);

///*do something*/
    int nbLabels = 256;
        GCoptimizationGridGraph gc(MXS,MYS,nbLabels);
  /* */
   
    for(y=0;y<MYS;y++)
    {
          for(x=0;x<MXS;x++)
          {
              int pixelIndex = y*MXS+x;
         
              for(int label=0;label<256;label++){
                   
                  gc.setDataCost(pixelIndex,label,(int)pow(label*255.0/nbLabels- image[pixelIndex],2));
               
              }

          }
    }
        
       /* Specifier le terme a priori */
    for ( int l1 = 0; l1 < nbLabels; l1++ )
            for (int l2 = 0; l2 < nbLabels; l2++ ){
                    float cost =  beta*abs(l1-l2);
                    gc.setSmoothCost(l1,l2,(int)cost);
            }     

  /* On lance l'algorithme de graph cut */
        gc.expansion(2);


 /* save results */
    int i=0;
    for( int y=0;y<MYS;y++){
        for( int x=0;x<MXS;x++,i++){
            c[i]=gc.whatLabel(i);
        }
    }
        
}
