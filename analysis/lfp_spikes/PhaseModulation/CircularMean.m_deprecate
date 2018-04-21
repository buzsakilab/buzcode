function [mu, Kappa, pval, Rmean, delta, sigma,confDw,confUp] = CircularMean(theta,varargin)

% [theta, Kappa, pval, Rmean, delta] = CircularMean(t) ciricualr mean, dispersion
% and significance of a circular datasets
% 
% INPUTS:
% t: a column vector of angles (in radians) or a n x p matrix, n: samples /
% p:variables
% OUTPUT:
% theta: the mean direction
% Kappa; the concentration factor to be used with Von Mises distribution
% pval: the signficance of the mean direction against an uniformity null
% Rmean: the mean resultant length
% delta: the sample ciricualr dispersion
% hypothesis (with a Rayleigh test) 
% see Fisher N.I. Analysis of Circular Data p. 30-35

% copyright (c) 2005 Francesco P. Battaglia, 2007 Adrien Peyrache
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html



% Default parameters
%B=200; %number of iterations in bootstrap
alpha = 0.05; %significancy level
confidence=0;


n = size(theta,1);
if n==0
	warning('No Data')
     mu = [];
    Kappa = [];
    pval = [];
    Rmean = [];
    delta = [];
    sigma = [];
    confDw = [];
    confUp = [];
    return
end

p = size(theta,2);
if p>1
    mu = [];
    Kappa = [];
    pval = [];
    Rmean = [];
    delta = [];
    sigma = [];
    confDw = [];
    confUp = [];
    for ii=1:p
        [mu0, Kappa0, pval0, Rmean0, delta0, sigma0,confDw0,confUp0] = CircularMean(theta(:,ii),varargin);
        mu = [mu;mu0];
        Kappa = [Kappa;Kappa0];
        pval = [pval;pval0];
        Rmean = [Rmean;Rmean0];
        delta = [delta;delta0];
        sigma = [sigma;sigma0];
        confDw = [confDw;confDw0];
        confUp = [confUp;confUp0];
    end
else


    while length(varargin)>1

        switch varargin{1}
            case 'alpha'
                if isnumeric(varargin{2})
                    alpha = varargin{2};
                    confidence=1;
                else
                    warning(['Argument must be numeric'])
                end
            otherwise
                warning(['Unknow option name, ' varargin{1}])
        end

        varargin = varargin(3:end);

    end

    S = sum(sin(theta));
    C = sum(cos(theta));

    mu = angleMean(C,S);
    Rmean = sqrt(S^2+C^2) / n;
    rho2 = sum(cos(2* (theta-mu))) / n;
    delta = (1 - rho2) / (2 * Rmean^2);
    sigma = sqrt(-2*log(Rmean));


    Z = n * Rmean^2;
    pval = exp(-Z) * (1 + (2*Z-Z^2)/(4*n) - (24*Z-132*Z^2+76*Z^3-9*Z^4)/(288*n^2));

    if Rmean <0.53 
        Kappa = 2*Rmean + Rmean^3 + 5*Rmean^5/6;
    elseif Rmean>=0.53 & Rmean<0.85
        Kappa = -0.4 + 1.39*Rmean + 0.43/(1-Rmean);
    else
        Kappa = 1/(Rmean^3 -4*Rmean^2 + 3*Rmean);
    end

    if n<15
        if Kappa<2
            Kappa = max(Kappa-2/(n*Kappa),0);
        else
            Kappa = (n-1)^3*Kappa/(n^3+n);
        end
    end


    confDw=0;
    confUp=0;

    if confidence

        x = cos(theta);
        y = sin(theta);
        z0 = mean([x y])';
        u0 = cov([x y]);
        v0 = sqrtm(u0);

        for i=1:B

            perm = floor(n*rand(n,1))+1;
            phi = theta(perm);
            x = cos(phi);
            y = sin(phi);
            zB = mean([x y])';
            uB = cov([x y]);
            vB = sqrtm(uB);
            wB = inv(vB);
            csB = z0 + v0*wB*(zB-z0);
            csB = csB/norm(csB);
            muB(i) = angleMean(csB(1),csB(2));

        end

        gamma = muB - mu;
        gamma = sort(gamma);
        l = floor(0.5*(B*alpha+1));
        m = B - l;
        confDw = mu + gamma(l+1);
        confUp = mu + gamma(m);

    end



end

end

function mu = angleMean(C,S)

%  	keyboard

	mu = atan2(S, C);
	if mu<0
		mu = mu+2*pi;
	end

end