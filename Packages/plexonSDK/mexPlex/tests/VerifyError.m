function  [results] = VerifyError(expr, results)
wasError = 0;
try 
    eval(expr);
catch ex
    wasError=1; 
end
if wasError == 0
    disp ([expr ' <-- expected error, but was no error'])
    count = 0;
    if(isfield(results, 'failedTests'))
        count = size(results.failedTests, 1);
    end
    results.failedTests{count+1,1} = expr;
    results.allTestsPassed = 0;
end
return