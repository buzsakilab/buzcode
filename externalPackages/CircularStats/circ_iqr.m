function p = circ_iqr(alpha)
% CIRC_IQR - estimates the interquartile range of circular data
%
% INPUT & OUTPUT
% -------------
% alpha - input data (in radians, not binned)
% p     - estimated iqr
%
% METHOD
% -------------
% For normal distributions IQR can be estimated from the 
% standard deviation of the sample, since 50% of the data
% lies within +/- 0.6745 sigma.
%     
% REQUIREMENTS
% -------------
% The functions uses the CircStat toolbox for circular data
% by Philipp Berens. Download from:
% http://www.mathworks.com/matlabcentral/fileexchange/10676
% 
% CREDITS
% -------------
% Vlad Atanasiu
% 2010.12.22
%     - creation

p = circ_std(alpha)/0.7413;
