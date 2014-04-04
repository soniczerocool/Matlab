
/*********************************************************************
 *
 *
 * MRF with data terms defined in a matrix.
 * 3D data
 *
 *
 * Keep in mind:
 * <> Use 0-based indexing as always in C or C++
 * <> Indexing is column-based as in Matlab (not row-based as in C)
 * <> Use linear indexing.  [x*dimy+y] instead of [x][y]
 *
 *
 ********************************************************************/
#include <matrix.h>
#include "C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\GCoptimization.h"
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
    mxArray *image_in_m,  *c_out_m , *matrix_in_m;
    const mwSize *dims , *size;
    double *image, *c, *matrix, beta;
    int MXS, MYS, MZS, numdims;
    int x,y,z;
    int Nb_points, nbLabels ;
    
    /* Check for proper number of input and output arguments */
    if (nrhs != 3) {
        mexErrMsgIdAndTxt( "MATLAB: markov_network_multilable:invalidNumInputs",
                "Wrong number of arguments\n\n\t Data usage : labelfield = demo2(image,mask,sigma).");
    }
    
    /* Check data type of input argument  */
    if (!(mxIsDouble(prhs[0])) || !(mxIsDouble(prhs[1]))  || !(mxIsDouble(prhs[2]))){
        mexErrMsgIdAndTxt( "MATLAB: markov_network_multilable:inputNotDouble",
                " Input arguments must be of type double.");
    }
    
    
//associate inputs
    image_in_m = mxDuplicateArray(prhs[0]);
    
    //d_in_m = mxDuplicateArray(prhs[2]);
//figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    MZS = (int)dims[0]; MYS = (int)dims[1];  MXS = (int)dims[2];
    
//associate outputs
    c_out_m = plhs[0] = mxCreateDoubleMatrix(MYS*MXS*MZS,1,mxREAL);
    // c_out_m = plhs[0] = mxCreateNumericArray(numdims, dims, mxINT16_CLASS, mxREAL);
    beta = mxGetScalar(prhs[2]);
//associate pointers
    image = mxGetPr(image_in_m);
    matrix_in_m = mxDuplicateArray(prhs[1]);
    matrix = mxGetPr(matrix_in_m);
    
    c = mxGetPr(c_out_m);
    
    size = mxGetDimensions(prhs[1]);
    Nb_points = (int)size[1],   nbLabels = (int)size[0];
    char buffer [33];
    
///*do something*/
    GCoptimizationGeneralGraph gc(MXS*MYS*MZS,nbLabels);
    
    int ii;
  for(z=0;z<MZS;z++)
    {  
  for(x=0;x<MXS;x++){
      for(y=0;y<MYS;y++)
    {
                unsigned int pixelIndex =z*(MXS*MYS)+ x*MYS+y;
                for(int label=0;label<nbLabels;label++)
                {
                    double datacost = matrix[pixelIndex * nbLabels + label];
                    gc.setDataCost(pixelIndex,label,datacost);
                }
                
            }
        }
    }
    /* Connect neighboring pixels */
for(z=1;z<MZS;z++){
   for(x=1;x<MXS;x++){
      for(y=1;y<MYS;y++){
               unsigned int pixelIndex1 = z*MXS*MYS+x*MYS+y;
                unsigned int pixelIndex2 = z*MXS*MYS+x*MYS+ y-1;
                unsigned int pixelIndex3 = z*MXS*MYS+ (x-1)*MYS+y;
                unsigned int pixelIndex4 = (z-1)*MXS*MYS+x*MYS+y; 
                gc.setNeighbors(pixelIndex1,pixelIndex2,beta);
                gc.setNeighbors(pixelIndex1,pixelIndex3,beta);
                gc.setNeighbors(pixelIndex1,pixelIndex4,beta);
                
            }
        }
    }
    
    
    /* go graph cut */
    gc.expansion(2);
    
    /* save results */
    int i=0;
    for( int z=0;z<MZS;z++) {
        for( int y=0;y<MYS;y++){
            for( int x=0;x<MXS;x++,i++){
                c[i]=gc.whatLabel(i);
            }
        }
    }
}
