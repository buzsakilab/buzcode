void mexmetypecheck(const mxArray *arr,int expected,char *errorMsg)
{
    int thetype = mxGetClassID(arr);
    char passed = 0;
    if(expected == mxINT8_CLASS || 
       expected == mxINT16_CLASS ||
       expected == mxINT32_CLASS ||
       expected == mxINT64_CLASS)
    {
        if(thetype == mxINT8_CLASS || 
           thetype == mxINT16_CLASS ||
           thetype == mxINT32_CLASS ||
           thetype == mxINT64_CLASS)
        {
            if(thetype <= expected)
                passed = 1;
        }
    }
    else if(expected == mxUINT8_CLASS || 
            expected == mxUINT16_CLASS ||
            expected == mxUINT32_CLASS ||
            expected == mxUINT64_CLASS)
    {
        if(thetype == mxUINT8_CLASS || 
           thetype == mxUINT16_CLASS ||
           thetype == mxUINT32_CLASS ||
           thetype == mxUINT64_CLASS)
        {
            if(thetype <= expected)
                passed = 1;
        }
    }
    else if(expected == mxSINGLE_CLASS || expected == mxDOUBLE_CLASS)
    {
        passed = 1;
        if(thetype == mxDOUBLE_CLASS && expected == mxSINGLE_CLASS)
        {
            passed = 0;
        }
    }

    if(!passed)
    {
        mexErrMsgTxt(errorMsg);
    }
}
