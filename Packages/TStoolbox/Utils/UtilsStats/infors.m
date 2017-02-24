% [I, It, sp] = infors(S, R, num_s, infors_opt)
% 
% MEX FILE
% 
% INPUTS:
%     S = a row or column vector containing the stimulus for each time bin
%     R = a num_trials x num_cells matrix containing the firing of each cell
%     during each time bin.
%     num_s = the number of different stimuli in S
%     infors_opt: a structure containing the options to be passed to infors:
%       seed: the seed for the random number generator
%       max_t: the max number of trials to be considered per each stimulus
%       min_t: the min number of trials to be considered per each stimulus
%       s0err: if non-zero, the stimulus with zero index is discarded
%       s_throw: throw away n stimuli to subsample
%       timewin: time window size in msec.
%       max_b: maximum number of spikes in a bin 
% 
% OUTPUTS:
%     I = a vector of mutual information for each cell in bits 
%     It = Information time derivative
%     sp = sparseness
% Treves, Panzeri & Skaggs 95/95
% mexified batta 1999
% status: BETA
% version 1