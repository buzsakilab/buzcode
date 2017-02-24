function SI = PlaceFields_StabilityIndex(fr1,fr2)

% USAGE:
%     SI = PlaceFields_StabilityIndex(fr1,fr2)
%     
% Computes stability index between two place maps.
% The index is definde as the average Hellinger distance between each bin, 
% assuming they are drawn from Poisson distributions.
%
% INPUTS:
%     fr1,2: place maps in conditions 1 and 2
% 
% OUTPUTS:
%     SI: stability index (between 0 (min) and 1 (max))
    
% Adrien Peyrache, 2015

coef = max(fr1,fr2);
mx = sum(coef(:));

H = coef/mx.*sqrt(1-exp(-0.5*(sqrt(fr1)-sqrt(fr2)).^2));

H = sum(H(:));

SI = 1 - H;