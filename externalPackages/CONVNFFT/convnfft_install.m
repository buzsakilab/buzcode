function convnfft_install
% function convnfft_install
% Installation by building the C-mex file needed for convnfft
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History
%  Original: 16/Sept/2009

arch=computer('arch');
mexopts = {'-O' '-v' ['-' arch]};

if ~verLessThan('MATLAB','9.4')
    R2018a_mexopts = {'-R2018a'};
else
    % 64-bit platform
    if ~isempty(strfind(computer(),'64'))
        mexopts(end+1) = {'-largeArrayDims'};
    end
    R2018a_mexopts = {};
end

% invoke MEX compilation tool
mex(mexopts{:},R2018a_mexopts{:},'inplaceprod.c');