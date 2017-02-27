function P=ar1spectrum(ar1,period)
% AR1 power spectrum... 
%
% power=ar1spectrum(ar1,period)
%
%
% (c) Aslak Grinsted 2002-2004
%


% -------------------------------------------------------------------------
%   Copyright (C) 2002-2004, Aslak Grinsted
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.


freq=1./period;
P=(1-ar1.^2)./(abs(1-ar1.*exp(-2*pi*i*freq))).^2; %http://www.madsci.org/posts/archives/may97/864012045.Eg.r.html
%fixed typo in numerical recipes