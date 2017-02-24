function [p,F] = CircularANOVA(angles,factors,method)

%CircularANOVA - One or two-way ANOVA on circular data.
%
%  Circular ANOVA can be performed using one of three methods with different
%  assumptions and limitations:
%
%    1) Watson & Williams (1956)
%         The data are drawn from von Mises distributions with equal concent-
%         ration parameters, the significance levels rely on X² approximations,
%         the extension to multi-way analysis is difficult to interpret because
%         interaction factors can be negative
%
%    2) Harrison & Kanji (1986)
%         The data are drawn from von Mises distributions with equal concent-
%         ration parameters, the significance levels rely on X² approximations,
%         the test can be extended to multi-way analysis (positive interaction
%         factors), the statistics are affected by both location and dispersion
%
%    3) Anderson & Wu (1995)
%         No distributional assumptions, the significance levels are determined
%         using a bootstrap approach, the test can be extended to multi-way
%         analysis
%
%  USAGE
%
%    [p,F] = CircularANOVA(angles,groups,method)
%    [p,F] = CircularANOVA(angles,[factor1 factor2],method)
%
%    angles         angles in radians
%    groups         group indices
%    factor1        levels for factor 1
%    factor2        levels for factor 2
%    method         method to compute ANOVA:
%                    'ww'     Watson & Williams (1956) (default)
%                    'l2'     Harrison & Kanji (1988), L2 distance
%                    'lr'     Anderson & Wu (1995), likelyhood ratio
%
%  SEE
%
%    See also CircularVariance, CircularConfidenceIntervals, Concentration,
%    ConcentrationTest.
%
%  REFERENCE
%
%    Anderson & Wu (1992) Decomposition of variation in directional data,
%    IIQP Report RR-92-10.

% Copyright (C) 2009-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

warning('This function is not fully implemented yet (not all methods are available).');

p = [];
F = [];
if isempty(angles), return; end

nIterations = 500;

% Check parameters
if nargin < 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help CircularANOVA">CircularANOVA</a>'' for details).');
end
if ~isdvector(angles),
	error('Incorrect angles (type ''help <a href="matlab:help CircularANOVA">CircularANOVA</a>'' for details).');
end
isradians(angles);

if nargin == 2,
	method = 'ww';
else
	method = lower(method);
	if ~isstring_FMAT(method,'ww','l2','lr'),
		error(['Unknown ''' method ''' method (type ''help <a href="matlab:help CircularANOVA">CircularANOVA</a>'' for details).']);
	end
end

A = exp(j*angles);
R = abs(mean(A));
N = length(A);
mu = angle(mean(A));

if isvector(factors),
	% One-way ANOVA
	groups = factors(:);
	groupIDs = unique(groups);
	if any(diff(groupIDs)~=1),
		error('Incorrect groups numbers (should be 1..N).');
	end
	switch(method),
		case 'ww',
			% Watson & Williams (1956)
			q = length(groupIDs);
			for i = 1:length(groupIDs),
			ok = groups == groupIDs(i);
				Ri(i) = abs(mean(A(ok)));
				Ni(i) = sum(ok);
				ki(i) = Concentration(angles(ok));
			end
			SSw = N - sum(Ni.*Ri); % ignore 2k factor
			SSb = sum(Ni.*Ri) - N*R; % ignore 2k factor
			F = (N-q)/(q-1)*SSb/SSw;
			k = sum(ki.*Ni)/N;
			if 2 < k && k < 10,
				% Improve X² approximation for moderate concentration (Stephens, 1969)
				F = F * (1+3/(8*k));
			end
			p = 1-fcdf(F,q-1,N-q);
		otherwise,
			error('One-way ANOVA: only the ''ww'' method is implemented.');
	end
else
	% Two-way ANOVA
	levels1 = unique(factor1);
	if any(diff(levels1)~=1),
		error('Incorrect levels for factor 1 (should be 1..N).');
	end
	levels2 = unique(factor2);
	if any(diff(levels2)~=1),
		error('Incorrect levels for factor 2 (should be 1..N).');
	end
	switch(method),
		case 'lr',
			% Anderson & Wu (1995), location-only
			% Check requirements: 2x2 balanced design
			if length(levels1)~=2,
				error('Two-way ANOVA: factors can only have 2 levels (2x2 design).');
			end
			if length(levels2)~=2,
				error('Two-way ANOVA: factors can only have 2 levels (2x2 design).');
			end
			m = Accumulate([factor1 factor2]);
			if any(m~=m(1)),
				error('Two-way ANOVA: all cells must have the same count (balanced design).');
			end
			m = m(1);
			% F values for observations
			[Fa,Fb,Fab] = lo(A,m,factor1,factor2);
			% F values for surrogate data
			for i = 1:nIterations,
				p = randperm(N);
				[sFa(i),sFb(i),sFab(i)] = lo(A(p),m,factor1,factor2);
			end
			p(1) = 1-npcdf(sFa,Fa);
			p(2) = 1-npcdf(sFb,Fb);
			p(3) = 1-npcdf(sFab,Fab);
			F = [Fa;Fb;Fab];
		otherwise,
			error('Two-way ANOVA: only the ''lr'' method is implemented.');
	end
end

%-----------------------------------------------------------------------------------------------------------------
%  Compute F values for 2-way ANOVA (location only method)
%-----------------------------------------------------------------------------------------------------------------

function [Fa,Fb,Fab] = lo(A,m,factor1,factor2)

% Mean angle
theta = angle(mean(A));

% Mean angles for each level of factor 1
for i = 1:length(levels1),
	theta_i(i) = angle(mean(A(factor1==i)));
end

% Mean angles for each level of factor 2
for i = 1:length(levels2),
	theta_j(i) = angle(mean(A(factor2==i)));
end

% Mean angles for each level of factor 1 and factor 2
for i1 = 1:length(levels1),
	for i2 = 1:length(levels2),
		theta_ij((i1-1)*length(levels2)+i2) = angle(mean(A(factor1==i1&factor2==i2)));
	end
end
t_ij = repmat(theta_ij,1,m);
t_ij = t_ij';
t_ij = t_ij(:);

% Interaction angles
theta_11 = mean(A(factor1==1&factor2==1));
theta_22 = mean(A(factor1==2&factor2==2));
theta_12 = mean(A(factor1==1&factor2==2));
theta_21 = mean(A(factor1==2&factor2==1));
phi_l(1) = angle(mean([theta_11 theta_22]));
phi_l(2) = angle(mean([theta_12 theta_21]));

% Sum of squares
SSa = 4*m*sum(1-cos(tehta_i-theta));
SSb = 4*m*sum(1-cos(tehta_j-theta));
SSab = 4*m*sum(1-cos(phi_l-theta));
SSr = 2*sum(1-cos(angles-t_ij));

% F values
Fa = SSa/SSr;
Fb = SSb/SSr;
Fab = SSab/SSr;
