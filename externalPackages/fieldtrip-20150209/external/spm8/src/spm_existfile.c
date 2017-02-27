/*
 * $Id: spm_existfile.c 8776 2013-11-14 09:04:48Z roboos $
 * Guillaume Flandin
 */
 
#define _FILE_OFFSET_BITS 64

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    bool status    = false;
    char *filename = NULL;
    FILE *fid      = NULL;
    
    if (nrhs != 1)
        mexErrMsgTxt("One input only required.");
    else
    {
        if (!mxIsChar(prhs[0]))
            mexErrMsgTxt("Input must be a string.");
        filename = mxArrayToString(prhs[0]);
        fid = fopen(filename,"r");
        if (fid != NULL)
        {
            status = true;
            fclose(fid);
        }
        mxFree(filename);
    }
        
    plhs[0] = mxCreateLogicalScalar(status);
}
