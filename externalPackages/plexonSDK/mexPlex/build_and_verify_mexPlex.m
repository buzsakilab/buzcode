function [res] = build_and_verify_mexPlex()
% this script will build and verify mexPlex
%
% OUTPUT:
%   returns 1 if all verification tests passed
%   returns 0 if some of the verification tests failed
% 
% if ReadingPLXandDDTfilesinMatlab-mexw.zip was inzipped to
%
% C:\PlexonMatlab
%
% this script should be run when Matlab's current directory is
%
% C:\PlexonMatlab\OfflineSDK\mexPlex
%
% If res is 1, all the tests passed.
%
% If you want to build mexPlex without verification, run this in Matlab:
%
% mex -output mexPlex -outdir .. PlexMethods.cpp
%
%
% Plexon made sure that mexPlex code compiles without warnings in Matlab R2009b 
% under Ubintu 10.04 using gcc 4.1 
% If you have trouble compiling mexPlex in Linux, please contact Matlab technical support.


% compile mexPlex
mex -output mexPlex -outdir .. PlexMethods.cpp

% make sure that <>\OfflineSDK is in Matlab path
cd ..
if length(findstr(path,pwd)) == 0
    addpath(pwd)
end
% go to tests
cd mexPlex
cd tests
% load the data extracted with mexPlex in Windows
load('mexPlexData1.dat', '-mat');

% verify that the newly built mexPlex generates the same data
t=evalc('res = verify_mexplex(data, pwd);');

if res == 0
    disp('VERIFICATION FAILED');
    t
else
    disp('VERIFICATION TESTS PASSED');
end
cd ..
