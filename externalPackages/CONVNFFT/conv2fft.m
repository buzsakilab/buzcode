function C = conv2fft(varargin)
% C = conv2fft(A, B)
% C = conv2fft(H1, H2, A)
%
%   C = CONV2FFT(A, B) performs the 2-D convolution of matrices A and B.
%   If [ma,na] = size(A), [mb,nb] = size(B), and [mc,nc] = size(C), then
%   mc = max([ma+mb-1,ma,mb]) and nc = max([na+nb-1,na,nb]).
%
%   C = CONV2FFT(H1, H2, A) first convolves each column of A with the vector
%   H1 and then convolves each row of the result with the vector H2.  If
%   n1 = length(H1), n2 = length(H2), and [mc,nc] = size(C) then
%   mc = max([ma+n1-1,ma,n1]) and nc = max([na+n2-1,na,n2]).
%   CONV2(H1, H2, A) is equivalent to CONV2FFT(H1(:)*H2(:).', A) up to
%   round-off.
%
%   C = CONV2FFT(..., SHAPE) returns a subsection of the 2-D
%   convolution with size specified by SHAPE:
%     'full'  - (default) returns the full 2-D convolution,
%     'same'  - returns the central part of the convolution
%               that is the same size as A.
%     'valid' - returns only those parts of the convolution
%               that are computed without the zero-padded edges.
%               size(C) = max([ma-max(0,mb-1),na-max(0,nb-1)],0).
%
%   See also CONV2, CONVN, CONVNFFT
%
%   Author: Bruno Luong <brunoluong@yahoo.com>
%   History:
%       Original: 21-April-2014

if length(varargin)>=3 && isnumeric(varargin{3})
    [H1, H2, A] = deal(varargin{1:3});
    varargin = varargin(4:end);
    C = convnfft(H1(:), A, varargin{:});
    C = convnfft(H2(:).', C, varargin{:});
else
    C = convnfft(varargin{:});
end

end % conv2fft
