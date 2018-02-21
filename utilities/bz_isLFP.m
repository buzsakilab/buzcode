function [isLFP] = bz_isLFP(lfp)
% USAGE
% [] = bz_isLFP(lfp)
% 
% INPUT
%       lfp   - struct with the following fields
%                   .lfp
%                   .timestamps
%                   .samplingRate
%                   .channels
%                   
%
% OUTPUT
%      logical true if struct meets lfp criteria, false if otherwise
%
% written by david tingley, 2017


if isfield(lfp,'data') && isfield(lfp,'timestamps') && isfield(lfp,'samplingRate') && isfield(lfp,'channels') % check that fields exist
     if ismatrix(lfp.data) && isvector(lfp.timestamps) && isnumeric(lfp.samplingRate) && isvector(lfp.channels)
         isLFP = true;
     else
         warning('one of the required fields for a lfp type is not formatted correctly')
         isLFP = false;
     end
else
    warning('one of the required fields for a lfp type does not exist')
    isLFP = false;
end
