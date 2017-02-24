function [TT, TTbad] = LoadTT(fn, offset, chunk_length)

% TT = LoadTT(fn)
%
% Loads an NSMA TT file.
%
% INPUTS:
%   fn -- .tt file
%   offset -- reading from this spike
%   chunk_length -- how many spikes to be read
%
% OUTPUTS:
%   TT is a tsd structure 
%      where data = nSpikes x nSamplesPerSpike x nTrodes
%
% uses mex file LoadTT0(fn) to do the main read.
%
% ADR 1998
% version L5.0
% status PROMOTED

if nargout == 1
  if nargin == 3
    [t, wv] = LoadTT1(fn, offset, chunk_length);
  elseif nargin == 1
    [t, wv] = LoadTT1(fn);
  end
  TT = tsd(t, wv);
elseif nargout == 2
    if nargin == 3
    [t, wv, tbad, wvbad] = LoadTT1(fn, offset, chunk_length);
  elseif nargin == 1
    [t, wv, tbad, wvbad] = LoadTT1(fn);
  end
  TT = tsd(t, wv);
  TTbad = tsd(tbad, wvbad);
end













